#include "FE_sweeper.hpp"



#include <dune/istl/bvector.hh> // BlockVector
#include <dune/istl/bcrsmatrix.hh> // BCRSMatrix

#include <dune/grid/yaspgrid.hh> // YaspGrid
#include <dune/functions/functionspacebases/pq1nodalbasis.hh> // PQ1NodalBasis
#include <dune/fufem/assemblers/dunefunctionsoperatorassembler.hh> // DuneFunctionsOperatorAssembler
#include <dune/fufem/assemblers/istlbackend.hh> // istlMatrixBackend
#include <dune/fufem/assemblers/localassemblers/massassembler.hh> //MassAssembler
#include <dune/fufem/assemblers/localassemblers/laplaceassembler.hh> //LaplaceAssembler
#include <dune/istl/paamg/graph.hh> // Dune::Amg::MatrixGraph
#include<dune/istl/matrixredistribute.hh> // Dune::RedistributeInformation
#include<dune/istl/preconditioners.hh> // Dune::BlockPreconditioner, Dune::SeqSSOR
#include<dune/istl/solvers.hh> // Dune::CGSolver, Dune::RestartedGMResSolver
#include<dune/istl/schwarz.hh> // Dune::OverlappingSchwarzScalarProduct, Dune::OverlappingSchwarzOperator


#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/istl/owneroverlapcopy.hh>// OwnerOverlapCopyCommunication

#include "mpi.h"


//this is a specified sweeper for the generalized Fisher equation
//it inherits from the implicit Finite Element sweeper

//you just need to add a few things

//the konstructer: mass- and stiffnessmatrix are already assembeld in the base class, but here we make a another discretisation vector w, because we use lumping for the nonlinearity

//exact: gives the exact solution of Fischers equation, we use it in the main to construct the initial value and to calculate the error after the simulation

//evaluate_rhs_impl: it gives back the right hand side of the discrtised ODE

//implicit_solve: it solves the resulting algebraic system which results from using sdc 

using namespace pfasst::examples::FE_sweeper;


namespace pfasst
{
  namespace examples
  {
    namespace fischer_example
    {
      template<
        class SweeperTrait,
        class BaseFunction,
        typename Enabled = void
      >
      class fischer_sweeper
        : public Heat_FE<SweeperTrait, BaseFunction, Enabled>{
            
        std::shared_ptr<VectorType>                     w; 
        double                                     	_nu{1.2};
        double                                     	_n{2.0};
        double                                      	_delta{1.0};
        double                                          _abs_newton_tol=1e-10; 

	MPI_Comm comm_x; 
            
	bool sequential = false;

        public:
            explicit fischer_sweeper<SweeperTrait, BaseFunction, Enabled>(std::shared_ptr<BaseFunction> basis, size_t nlevel, std::shared_ptr<GridType> grid)
                                    : Heat_FE<SweeperTrait, BaseFunction, Enabled>(basis, nlevel, grid){
        
                
                w = std::make_shared<VectorType>(this->M_dune.M());
        
                for(int j=0; j<this->M_dune.M(); j++){(*w)[j]=0;}

                for(int i=0; i<this->M_dune.M(); i++){
                    for(int j=0; j<this->M_dune.M(); j++){
                        if(this->M_dune.exists(i,j))
                            (*w)[i][0]= ((double) (*w)[i][0]) + ((double) this->M_dune[i][j][0][0]);
                    }
                }


          }
            
            
            
        
          fischer_sweeper(const fischer_sweeper<SweeperTrait, BaseFunction, Enabled>& other) = default;
          fischer_sweeper(fischer_sweeper<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          virtual ~fischer_sweeper() = default;

	  void set_comm(MPI_Comm comm){comm_x=comm;}          
  
          shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t)      {
            auto result = this->get_encap_factory().create();
            const auto dim = 1; //SweeperTrait::DIM;
            double n  = this-> _n;
            double l0 = this-> _nu;
            double l1 = l0/2. *(pow((1+n/2.), 1/2.) + pow((1+ n/2.), -1/2.) );
            double d = l1 - pow(pow(l1,2) - pow(l0,2), 1/2.);
            auto exact_solution = [l0, l1, n, d, t](const Dune::FieldVector<double,dim>&x){ 
                return pow((1 + (pow(2, n/2.)-1 )* exp(-(n/2.)*d*(x+2*l1*t)) ), -2./n);
            };  
            auto N_x = [t](const Dune::FieldVector<double,dim>&x){
                return x;
            };
            Dune::BlockVector<Dune::FieldVector<double,dim>> x_node;
            interpolate(*this->basis, x_node, N_x);
            interpolate(*this->basis, result->data(), exact_solution);
            return result;
          }
          
          shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& t,const shared_ptr<typename SweeperTrait::encap_t> u) {
	
            ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << t);
            auto result = this->get_encap_factory().create();
            double nu =this->_nu;
            for (int i=0; i<u->get_data().size(); ++i)
                result->data()[i]= -pow(u->get_data()[i], this->_n+1) * (*w)[i];	
            this->M_dune.umv(u->get_data(), result->data());
            result->data()*=nu*nu;
            this->A_dune.umv(u->get_data(), result->data());
        
            return result;
        }
        
        void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& t,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs){
	
         ML_CVLOG(4, this->get_logger_id(), "IMPLICIT spatial SOLVE at t=" << t << " with dt=" << dt);

         auto residuum = this->get_encap_factory().create();
	 Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 , newton_rhs_seq;
         newton_rhs.resize(rhs->get_data().size());
         newton_rhs2.resize(rhs->get_data().size());
         newton_rhs_seq.resize(rhs->get_data().size());
    
			for(int k=0; k< u->data().size(); k++){
		//std::cout << u->data()[k] << exact(0)->data()[k] << std::endl;
    		}

         u->zero();
	 for (int i=0; i< 200 ;i++){
            Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
            evaluate_f(f, u, dt, rhs);
            evaluate_df(df, u, dt);
            df.mv(u->data(), newton_rhs);
            newton_rhs -= f->data();
	    newton_rhs_seq = newton_rhs;	
            newton_rhs2 = newton_rhs;
		int rank;
		MPI_Comm_rank(comm_x, &rank );

	    auto u_seq = this->get_encap_factory().create();	

	    //if (sequential){	          
            	/*Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df);
	  
            	Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df,1.0);
	  
            	Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-16, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
          
          
            	Dune::InverseOperatorResult statisticss ;
            	cg.apply(u_seq->data(), newton_rhs_seq , statisticss ); */ //rhs ist nicht constant!!!!!!!!!
	    //}else{




		using DuneVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
		using DuneMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;

		MPI_Comm dune_comm = comm_x; //MPI_COMM_WORLD;
		
		using DuneCommunication = Dune::OwnerOverlapCopyCommunication<std::size_t>;
		DuneCommunication DuneComm(dune_comm);
		DuneCommunication *comm_redist;
		
		DuneMatrixType parallel_A;
		using DuneMatrixGraph = Dune::Amg::MatrixGraph<DuneMatrixType>;
		Dune::RedistributeInformation<DuneCommunication> dune_rinfo;

		bool hasDofs = Dune::graphRepartition(DuneMatrixGraph(df), 
										  DuneComm,
										  static_cast<int>(2),
										  comm_redist,
										  dune_rinfo.getInterface(),
										  true);

		dune_rinfo.setSetup();
		redistributeMatrix(df, parallel_A, DuneComm, *comm_redist, dune_rinfo);

	int size = u->data().size();
	/* displayinng the distributed part of matrix on each process*/
	/*for(int K=0; K<4; ++K){
	if(rank == K)
	{
		std::cout<<"Now Displaying the matrix part from process "<< rank <<std::endl;
		const int size = parallel_A.M();
		for(int i=0; i<size; ++i)	
		{
			for(int j=0; j<size; ++j)
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
	}*/

		comm_redist->remoteIndices().rebuild<false>();
		using Seq_Preconditioner = Dune::SeqSSOR<DuneMatrixType, DuneVectorType, DuneVectorType>;
				
		using Par_Preconditioner = Dune::BlockPreconditioner<DuneVectorType, DuneVectorType, DuneCommunication, Seq_Preconditioner>;
		using Par_ScalarProduct = Dune::OverlappingSchwarzScalarProduct<DuneVectorType, DuneCommunication>;
		using Par_LinearOperator = Dune::OverlappingSchwarzOperator<DuneMatrixType, DuneVectorType, DuneVectorType, DuneCommunication>;
	
		Par_ScalarProduct parallel_sp(*comm_redist);
		Par_LinearOperator parallel_linearoperator(parallel_A, *comm_redist);
		Seq_Preconditioner seq_precon(parallel_A, 1, 1.0); 
		Par_Preconditioner parallel_precon(seq_precon, *comm_redist);

		Dune::InverseOperatorResult statistics;

	
		const double residual_tol = 1e-16;
		const int num_restart = 10;		
		const int max_iter = 200;
		
		const int verbosity = rank ==0? 0:0;

		Dune::RestartedGMResSolver<DuneVectorType> GMRES(parallel_linearoperator, parallel_sp,
														 parallel_precon,
														 residual_tol,
														 num_restart,
														 max_iter,
														 verbosity);
		
		DuneVectorType parallel_b(parallel_A.N());
		DuneVectorType parallel_x(parallel_A.M());
		dune_rinfo.redistribute(newton_rhs, parallel_b);

		GMRES.apply(parallel_x, parallel_b, statistics);

		dune_rinfo.redistributeBackward(u->data(), parallel_x);
	
		if(rank==0) 
			for(int i=0; i< u->data().size(); i++)
				std::cout << " unterschied " << u->data()[i] << " " << u_seq->data()[i] << std::endl;
  
		//std::cout << rank << std::endl;
		MPI_Bcast(&(u->data()[0][0]), u->data().size(), MPI_DOUBLE, 0, comm_x);



	    //}	


	    //u_seq->data()  -= u->data();
	    //double n = u_seq->norm0();		
		
	    //u->data() = u_seq->data();	
            evaluate_f(f, u, dt, rhs);
          
            //std::cout << rank << " unterschied sequentiell  " << n << std::endl;  
	    //if(rank==0) std::exit(0);	
            if(f->norm0()<1e-10){   std::cout << "genauigkeit erreicht " << i << std::endl;      break;} //  std::exit(0); std::cout << "genauigkeit erreicht " << i << std::endl;
          
            df.mv(u->data(), residuum->data());
            residuum->data() -= newton_rhs2;
            //std::cout << "residuums norm " << residuum->norm0() << std::endl;

	}

	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }
        //evaluate_rhs_impl(0, u);
	//std::exit(0);
        this->_num_impl_solves++;
      }
      
        private:
            
      void evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
            const shared_ptr<typename SweeperTrait::encap_t> u,
            const typename SweeperTrait::time_t& dt,
            const shared_ptr<typename SweeperTrait::encap_t> rhs){

          double _nu=this->_nu;
          double _n=this->_n;

          f->zero();
          Dune::BlockVector<Dune::FieldVector<double,1> > fneu;
          fneu.resize(u->get_data().size());
          for (int i=0; i<u->get_data().size(); ++i)
          {
            f->data()[i]= pow(u->get_data()[i], _n+1) * (*w)[i];	
          }
          this->M_dune.mmv(u->get_data(), f->data());

          f->data() *= (_nu*_nu);


          this->A_dune.mmv(u->get_data(),f->data());
          f->data() *= dt;
          this->M_dune.umv(u->get_data(),f->data());
          f->data() -=rhs->get_data();
      }
						
      void evaluate_df(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u,
						 const typename SweeperTrait::time_t& dt
 						){
            double _nu=this->_nu;
            double _n=this->_n;
            for (int i=0; i<df.N(); ++i)
            {
                df[i][i]= (_nu*_nu)*(_n+1) * pow(u->get_data()[i], _n) * ((double) (*w)[i]);	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
      } 
          
    };   
    }
  }
}
