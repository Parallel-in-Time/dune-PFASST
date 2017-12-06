#include "FE_sweeper.hpp"



#include "assemble.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>

#include <iostream>


#include <dune/istl/bvector.hh> // BlockVector
#include <dune/istl/bcrsmatrix.hh> // BCRSMatrix

/*#include <dune/grid/yaspgrid.hh> // YaspGrid
#include <dune/functions/functionspacebases/pq1nodalbasis.hh> // PQ1NodalBasis
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh> // DuneFunctionsOperatorAssembler
#include <dune/fufem/assemblers/istlbackend.hh> // istlMatrixBackend
#include <dune/fufem/assemblers/localassemblers/massassembler.hh> //MassAssembler
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh> //LaplaceAssembler*/
#include <dune/istl/paamg/graph.hh> // Dune::Amg::MatrixGraph
#include <dune/istl/matrixredistribute.hh> // Dune::RedistributeInformation
#include <dune/istl/preconditioners.hh> // Dune::BlockPreconditioner, Dune::SeqSSOR
#include <dune/istl/solvers.hh> // Dune::CGSolver, Dune::RestartedGMResSolver
#include <dune/istl/schwarz.hh> // Dune::OverlappingSchwarzScalarProduct, Dune::OverlappingSchwarzOperator


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication


#include "mpi.h"




using std::shared_ptr;
using std::vector;








namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::init_opts()
      {
        /*config::options::add_option<size_t>("Heat FE", "num_dofs",
                                            "number spatial degrees of freedom per dimension on fine level");
        config::options::add_option<size_t>("Heat FE", "coarse_factor",
                                            "coarsening factor");
        config::options::add_option<spatial_t>("Heat FE", "nu",
                                               "thermal diffusivity");*/
      }


        template<class SweeperTrait, typename Enabled>
      Heat_FE<SweeperTrait, Enabled>::Heat_FE( std::shared_ptr<fe_manager> FinEl, size_t nlevel) // std::shared_ptr<fe_manager> FinEl,
        :   IMEX<SweeperTrait, Enabled>()

      {

	this->FinEl = FinEl;
	basis = FinEl->get_basis(nlevel);
	
		int rank, num;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank );    
	if (rank==0) assembleProblem(basis, this->A_dune, this->M_dune);

        const auto bs = basis->size();
        std::cout << "Finite Element basis consists of " <<  basis->size() << " elements " << std::endl;

        this->encap_factory()->set_size(bs);


      }




      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::set_options()
      {


        IMEX<SweeperTrait, Enabled>::set_options();

        this->_nu = config::get_value<spatial_t>("nu", this->_nu);

        std::cout << "thermal diffusivity is set to "   << this->_nu << std::endl;

        int num_nodes = this->get_quadrature()->get_num_nodes();



      }

        template<class SweeperTrait, typename Enabled>
        template<typename Basis>
      void
      Heat_FE<SweeperTrait, Enabled>::assemble(Basis &basis){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::assemble(){

        assembleProblem(basis, this->A_dune, this->M_dune);


      };




      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_t& t)
      {
        auto result = this->get_encap_factory().create();

        const auto dim = DIMENSION;
        spatial_t nu = this-> _nu;

        auto exact_solution = [t, nu, dim](const InVectorType &x){
            double solution=1.0;
            for(int i=0; i<dim; i++){solution *= std::sin(PI * x[i]);}
            return solution * std::exp(-t * dim * PI*PI * nu);
        };

        auto N_x = [t](const InVectorType &x){
            
            return x;

        };

        VectorType x_node;
        interpolate(*basis, x_node, N_x);

        interpolate(*basis, result->data(), exact_solution);
       

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        ML_CLOG(INFO, this->get_logger_id(), "number function evaluations:");
        //ML_CLOG(INFO, this->get_logger_id(), "  expl:        " << this->_num_expl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl:        " << this->_num_impl_f_evals);
        ML_CLOG(INFO, this->get_logger_id(), "  impl solves: " << this->_num_impl_solves);

        //this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged(const bool pre_check)
      {
        const bool converged = IMEX<SweeperTrait, Enabled>::converged(pre_check);

        if (!pre_check) {
          assert(this->get_status() != nullptr);
          const typename traits::time_t t = this->get_status()->get_time();
          const typename traits::time_t dt = this->get_status()->get_dt();

          assert(this->get_quadrature() != nullptr);
          auto nodes = this->get_quadrature()->get_nodes();
          const auto num_nodes = this->get_quadrature()->get_num_nodes();
          nodes.insert(nodes.begin(), typename traits::time_t(0.0));

          ML_CVLOG(1, this->get_logger_id(),
                   "Observables after "
                   << ((this->get_status()->get_iteration() == 0)
                          ? std::string("prediction")
                          : std::string("iteration ") + std::to_string(this->get_status()->get_iteration())));
          for (size_t m = 0; m < num_nodes; ++m) {
            ML_CVLOG(1, this->get_logger_id(),
                     "  t["<<m<<"]=" <<  (t + dt * nodes[m])
                     << "      |abs residual| = " <<  this->_abs_res_norms[m]
                     << "      |rel residual| = " <<  this->_rel_res_norms[m]
//                      << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[m])
//                      << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[m])
                    );
          }
          ML_CLOG(INFO, this->get_logger_id(),
                  "  t["<<num_nodes<<"]=" <<  (t + dt * nodes[num_nodes])
                  << "      |abs residual| = " << this->_abs_res_norms[num_nodes]
                  << "      |rel residual| = " << this->_rel_res_norms[num_nodes]
//                   << "      |abs error| = " << LOG_FLOAT << encap::norm0(error[num_nodes])
//                   << "      |rel error| = " << LOG_FLOAT << encap::norm0(rel_error[num_nodes])
                 );
        }
        return converged;
      }

      template<class SweeperTrait, typename Enabled>
      bool
      Heat_FE<SweeperTrait, Enabled>::converged()
      {
        return this->converged(false);
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat_FE<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory().size();
      }



      template<class SweeperTrait, typename Enabled>
      shared_ptr<GridType>
      Heat_FE<SweeperTrait, Enabled>::get_grid() const
      {
        return grid;
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_error(const typename SweeperTrait::time_t& t)
      {
        ML_CVLOG(4, this->get_logger_id(), "computing error");

        assert(this->get_status() != nullptr);
        const typename traits::time_t dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), typename traits::time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> error;
        error.resize(num_nodes + 1);
        std::generate(error.begin(), error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const typename traits::time_t ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), this->get_states()[m]);
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_t>>
      Heat_FE<SweeperTrait, Enabled>::compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                                            const typename SweeperTrait::time_t& t)
      {
        UNUSED(t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_t(0.0));

        vector<shared_ptr<typename traits::encap_t>> rel_error;
        rel_error.resize(error.size());
        std::generate(rel_error.begin(), rel_error.end(),
                 std::bind(&traits::encap_t::factory_t::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          rel_error[m]->scaled_add(1.0 / this->get_states()[m]->norm0(), error[m]);
        }

        return rel_error;
      }

      /*template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {
        UNUSED(u);
        ML_CVLOG(4, this->get_logger_id(),  "evaluating EXPLICIT part at t=" << t);
        auto result = this->get_encap_factory().create();
        result->zero();
        this->_num_expl_f_evals++;
        return result;
      }*/

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_t>
      Heat_FE<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                                                       const shared_ptr<typename SweeperTrait::encap_t> u)
      {



        ML_CVLOG(4, this->get_logger_id(), "evaluating IMPLICIT part at t=" << t);


        auto result = this->get_encap_factory().create();

        double nu =this->_nu;

        this->A_dune.mmv(u->get_data(), result->data());


        result->data() *= nu;
        //std::cout << "f_impl mit evaluate " << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          //std::cout << "f u " << result->data()[i] << std::endl;
        }

        
        return result;
        


      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat_FE<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs)
      {




	MPI_Comm comm_x=MPI_COMM_WORLD; 
 	//int rank, num;
	int rank, num;
	MPI_Comm_rank(comm_x, &rank );
	MPI_Comm_size(comm_x, &num );
	MPI_Barrier(comm_x);
	//std::cout << "test "<< rank << " von " << num << std::endl; //<< u_seq[i] << std::endl;
        Dune::BlockVector<Dune::FieldVector<double,1> > M_rhs_dune ;
        Dune::BlockVector<Dune::FieldVector<double,1> > u_seq;
        M_rhs_dune.resize(rhs->get_data().size());
	
	
        M_rhs_dune = rhs->get_data(); 

       
	
        Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune = 	Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->A_dune);
        M_dtA_dune *= (dt * this->_nu);
        M_dtA_dune += this->M_dune;


	
        /*if(rank==0){Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(M_dtA_dune);

        Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(M_dtA_dune,1.0);

        Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-10, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver

        Dune::InverseOperatorResult statisticss ;

        cg.apply(u_seq, M_rhs_dune , statisticss ); }*/ //rhs ist nicht constant!!!!!!!!!

    auto dune_comm = Dune::MPIHelper::getCollectiveCommunication();
		//MPI_Comm dune_comm = comm_x;//MPI_COMM_WORLD;
		int size = u->data().size();
/*if(rank==0)
for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
      		{
          		if(M_dtA_dune.exists(i,j))
          		{
	        		std::cout<<M_dtA_dune[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}

		

	MPI_Barrier(dune_comm);*/

	//start parallel


		
		using DuneVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
		using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
		using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>;
		DuneCommunication DuneComm(dune_comm);
		DuneCommunication *comm_redist;
		
		DuneMatrixType parallel_A;
		using DuneMatrixGraph = Dune::Amg::MatrixGraph<DuneMatrixType>;
		Dune::RedistributeInformation<DuneCommunication> dune_rinfo;

		bool hasDofs = Dune::graphRepartition(DuneMatrixGraph(M_dtA_dune), 
										  DuneComm,
										  static_cast<int>(num),
										  comm_redist,
										  dune_rinfo.getInterface(),
										  false);

 


		dune_rinfo.setSetup();
		redistributeMatrix(M_dtA_dune, parallel_A, DuneComm, *comm_redist, dune_rinfo);





		comm_redist->remoteIndices().rebuild<false>();

	/* displayinng the distributed part of matrix on each process*/
	/*for(int K=0; K<4; ++K){
	if(rank == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<< rank <<std::endl;
		const int si = parallel_A.N();
		for(int i=0; i<si; ++i)	
		{
			for(int j=0; j<si; ++j)
      		{
          		if(parallel_A.exists(i,j))
          		{
	        		std::cout<<parallel_A[i][j] << " ";
          		}
				else
				{
					std::cout<< 0 << " ";
				}
			}
			std::cout<<std::endl;
		}
	}	
	MPI_Barrier(dune_comm);
	}*/ //std::exit(0);


		using Seq_Preconditioner = Dune::SeqSSOR<DuneMatrixType, DuneVectorType, DuneVectorType>;
				
		using Par_Preconditioner = Dune::BlockPreconditioner<DuneVectorType, DuneVectorType, DuneCommunication, Seq_Preconditioner>;
		using Par_ScalarProduct = Dune::OverlappingSchwarzScalarProduct<DuneVectorType, DuneCommunication>;
		using Par_LinearOperator = Dune::OverlappingSchwarzOperator<DuneMatrixType, DuneVectorType, DuneVectorType, DuneCommunication>;
	
		Par_ScalarProduct parallel_sp(*comm_redist);
		Par_LinearOperator parallel_linearoperator(parallel_A, *comm_redist);
		Seq_Preconditioner seq_precon(parallel_A, 1, 1.0); 
		Par_Preconditioner parallel_precon(seq_precon, *comm_redist);

		Dune::InverseOperatorResult statistics;

		//std::cout << "test 1 "<< std::endl; //<< u_seq[i] << std::endl;
		const double residual_tol = 1e-10;
		const int num_restart = 200;		
		const int max_iter = 200;
		
		const int verbosity = 0;//rank ==0? 0:0;

		Dune::RestartedGMResSolver<DuneVectorType> GMRES(parallel_linearoperator, parallel_sp,
														 parallel_precon,
														 residual_tol,
														 num_restart,
														 max_iter,
														 verbosity);
		

		/*Dune::CGSolver<VectorType> GMRES(parallel_linearoperator, parallel_sp,
														 parallel_precon,
														 residual_tol,
														 max_iter, verbosity);*/



		DuneVectorType parallel_b(parallel_A.N());
		DuneVectorType parallel_x(parallel_A.M());
		dune_rinfo.redistribute(M_rhs_dune, parallel_b);

          //std::cout << "groesse "<< M_dtA_dune.N() << " " << parallel_A.N() << " " << parallel_A.M() <<" " << M_dtA_dune.M()<<  std::endl;
          //std::cout << "das ist jetzt das rhs " <<  std::endl;
	
	/*MPI_Comm_rank(MPI_COMM_WORLD, &rank );
	MPI_Barrier(MPI_COMM_WORLD);  	
	if(rank==0) for (size_t i = 0; i < parallel_b.size(); i++) {
          std::cout << "result " << parallel_b[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);
	if(rank==1) for (size_t i = 0; i < parallel_b.size(); i++) {
          std::cout << "result " << parallel_b[i] << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD); */


        /*std::cout << "das ist jetzt die matrix " <<  std::endl;	
	if(rank==0) for (size_t i = 0; i < parallel_A.M(); i++) {
			for (size_t j = 0; j < parallel_A.N(); j++){
          			if (parallel_A.exists(i,j)) {std::cout <<  parallel_A[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        }MPI_Barrier(MPI_COMM_WORLD);	if(rank==1)for (size_t i = 0; i < parallel_A.M(); i++) {
			for (size_t j = 0; j < M_dtA_dune.N(); j++){
          			if (parallel_A.exists(i,j)) {std::cout << parallel_A[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        } 
	MPI_Barrier(MPI_COMM_WORLD); 
	if(rank==2)for (size_t i = 0; i < parallel_A.M(); i++) {
			for (size_t j = 0; j < parallel_A.N(); j++){
          			if (parallel_A.exists(i,j)) {std::cout <<  parallel_A[i][j] << " ";}else{std::cout  << 0 << " ";} }std::cout << std::endl;
        } 
	MPI_Barrier(MPI_COMM_WORLD); */
	//std::exit(0);

	//std::exit(0);


//double starttime, endtime;
//       starttime = MPI_Wtime();

		GMRES.apply(parallel_x, parallel_b, statistics);

//endtime   = MPI_Wtime();
//       printf("That took %f seconds\n",endtime-starttime);
//std::exit(0);

		dune_rinfo.redistributeBackward(u->data(), parallel_x);
	
		/*if(rank==0) 
			for(int i=0; i< u->data().size(); i++)
				std::cout << " unterschied " << u->data()[i] << " "<< std::endl; //<< u_seq[i] << std::endl;*/
  
		//std::cout << rank << std::endl;
		MPI_Bcast(&(u->data()[0][0]), u->data().size(), MPI_DOUBLE, 0, comm_x);
	//end parallel









	


	
	
	
        ML_CVLOG(4, this->get_logger_id(),
                 "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);


	
        Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
        this->M_dune.mv(u->get_data(), M_u);
	

        //std::cout << "f_impl mit impl_solve" << std::endl;
        for (size_t i = 0; i < u->get_data().size(); i++) {
          //f->data()[i] = (u->data()[i] - rhs->data()[i]) / (dt);
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
          //std::cout << "f u " << f->data()[i] << std::endl;
        }

        this->_num_impl_solves++;
        //if (this->_num_impl_solves==5) std::exit(0);





      }
    }  // ::pfasst::examples::heat1
  }  // ::pfasst::examples
} // ::pfasst
