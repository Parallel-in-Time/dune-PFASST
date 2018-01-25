#include "sp.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"

#include <dune/common/identitymatrix.hh>

namespace pfasst
{
  namespace contrib
  {
  
    
 

    
    
    
    template<class TransferTraits>
    void
    SpectralTransfer<
            TransferTraits>
	::create(std::shared_ptr<fe_manager> FinEl, GFS& gfs0, GFS& gfs1)
    {
	
	      std::shared_ptr<std::vector<MatrixType*>> vecvec(FinEl->get_transfer());
	      std::shared_ptr<std::vector<MatrixType*>> vecvec_t(FinEl->get_transfer_t());

	      set_matrix(*vecvec->at(0), *vecvec_t->at(0), gfs0, gfs1);
	    
    }    
    
    
    
    

    template<class TransferTraits>
    void
    SpectralTransfer<
            TransferTraits>::set_matrix(Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> interpolate, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> restrict, GFS& gfs0, GFS& gfs1)
    {




	    interpolate_matrix = interpolate;
	    
	    restrict_matrix   = restrict;

	    inject_matrix   = restrict;

	    for (int i=0; i< restrict_matrix.N(); i++){
	      for (int j=0; j< restrict_matrix.M(); j++){
		if(inject_matrix.exists(i,j)){	
		  if (inject_matrix[i][j]==0.5 ) inject_matrix[i][j]=0;
		}

	      }
	    }    //injection fehlt!!!!!!!!!!!!!!!!

	//NOO novlpOperator(gfs_0, (*transferMatrix)[0][0]);

	//novlpOperator.apply(v0, v1);
		int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );






		/*if(my_rank==1)
		for(int i=0; i<interpolate_matrix.N(); i++){
		for(int j=0; j<interpolate_matrix.M(); j++){ 
			if (interpolate_matrix.exists(i,j)) {
				std::cout << interpolate_matrix[i][j]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl; 
		}


		if(my_rank==1)
		for(int i=0; i<restrict_matrix.N(); i++){
		for(int j=0; j<restrict_matrix.M(); j++){ 
			if (restrict_matrix.exists(i,j)) {
				std::cout << restrict_matrix[i][j]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl; 
		}

		if(my_rank==1)
		for(int i=0; i<inject_matrix.N(); i++){
		for(int j=0; j<inject_matrix.M(); j++){ 
			if (inject_matrix.exists(i,j)) {
				std::cout << inject_matrix[i][j]<< " ";
			}else { 
				std::cout << 0 << " ";
			} 
		}
		std::cout << std::endl; 
		}

	std::fill(n.begin(), n.end(), 40);
        gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); 

        gridp->globalRefine(1);

        fem = std::make_shared<FEM>((gridp->levelGridView(0))); 

  	gfs  = std::make_shared<GFS>(gridp->levelGridView(0),*fem); 


	prolong = std::make_shared<NOO>(*gfs, interpolate_matrix);
	rest = std::make_shared<NOO>(*gfs, restrict_matrix);
	inj =  std::make_shared<NOO>(*gfs, inject_matrix);
        std::cout <<  "+++ interpolate " <<  std::endl;



	using Zl = Dune::PDELab::Backend::Vector<GFS,double>;
	Zl v0(*gfs);

	//prolong = std::make_shared<NOO>(*gfs, interpolate_matrix);

	std::cout << "prolong gesetzt " << std::endl;


	/*Dune::TransposedMatMultMatResult< Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>>, Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> > m;
	Dune::IdentityMatrix< Dune::FieldMatrix<double, 1, 1>, 1 >  id;
	typedef Dune::FieldMatrix<double,1,1> M;
	Dune::BCRSMatrix<M> B(interpolate_matrix.N(),interpolate_matrix.M(), Dune::BCRSMatrix<M>::random);
	for (int i=0; i<interpolate_matrix.M();i++) B.setrowsize(i,1);
	B.endrowsizes();
	for (int i=0; i<interpolate_matrix.M();i++) B.addindex(i,i);;
	B.endindices();

        Dune::matMultTransposeMat(t_interpolate_matrix,  interpolate_matrix, B, true);*/


    }


    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_t> coarse,
                                          shared_ptr<typename TransferTraits::fine_encap_t> fine)
    {
      ML_CVLOG(1, "TRANS", "interpolate data");

      const size_t coarse_ndofs = coarse->get_data().N();
      const size_t fine_ndofs = fine->get_data().N();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        fine->data() = coarse->get_data();

      } else {
	/*const int dim=1;
	const int degree =1;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<double,dim> h = {1};	      
	std::array<int,dim> n;
	std::fill(n.begin(), n.end(), 200);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(h,n));

        gridp->refineOptions(false); 
        gridp->globalRefine(1);
	typedef Grid::LevelGridView GV;
        GV gv1=gridp->levelGridView(1);
        GV gv0=gridp->levelGridView(0);
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
        FEM fem0(gv0);
        FEM fem1(gv1);
  	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON;
  	typedef Dune::PDELab::istl::VectorBackend<> VBE;
  	typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  	GFS gfs0(gv0,fem0);*/

	//prolong = std::make_shared<NOO>(gfs0, interpolate_matrix);


	//Dune::BCRSMatrix <Dune::FieldMatrix<double, 1, 1>> test(interpolate_matrix.begin(), interpolate_matrix.end());
	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	MPI_Barrier(MPI_COMM_WORLD); 
	if (my_rank==0) Dune::PDELab::Backend::native(fine->data())[0][0] = 0;
	MPI_Barrier(MPI_COMM_WORLD); 
	if (my_rank==num_pro-1) Dune::PDELab::Backend::native(fine->data())[Dune::PDELab::Backend::native(fine->data()).N()-1][Dune::PDELab::Backend::native(fine->data()).N()-1] = 0;
	if(my_rank==0) for(auto r =coarse->data().begin(); r !=coarse->data().end(); ++r){std::cout << "coarse " << *r <<std::endl;} //std::exit(0);
 	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for(auto r =coarse->data().begin(); r !=coarse->data().end(); ++r){std::cout << "coarse " << *r <<std::endl;} //std::exit(0);
        //std::cout <<  "+++ vor prolong apply " <<  std::endl;
	//prolong->apply(coarse->data(), fine->data());
        //std::cout <<  "+++ nach prolong apply " <<  std::endl;
	if(my_rank==0){for(int i =0; i<  Dune::PDELab::Backend::native(fine->data()).N()-1; i++){
		if(i%2==0){
			Dune::PDELab::Backend::native(fine->data())[i] = Dune::PDELab::Backend::native(coarse->data())[i/2] ;
		}else{
			Dune::PDELab::Backend::native(fine->data())[i] = (Dune::PDELab::Backend::native(coarse->data())[i/2] + Dune::PDELab::Backend::native(coarse->data())[i/2+1] )*0.5;
		}
	}}else{for(int i =1; i<  Dune::PDELab::Backend::native(fine->data()).N(); i++){
		if(i%2==1){
			Dune::PDELab::Backend::native(fine->data())[i] = Dune::PDELab::Backend::native(coarse->data())[(i+1)/2] ;
		}else{
			Dune::PDELab::Backend::native(fine->data())[i] = (Dune::PDELab::Backend::native(coarse->data())[i/2] + Dune::PDELab::Backend::native(coarse->data())[i/2+1] )*0.5;
		}
	}}

        //interpolate_matrix.mv(Dune::PDELab::Backend::native(coarse->data()), Dune::PDELab::Backend::native(fine->data()));
	if(my_rank==0) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine i" << *r <<std::endl;}	MPI_Barrier(MPI_COMM_WORLD);  std::cout << "proc " << std::endl; MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==1) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine i" << *r <<std::endl;}// std::exit(0);
        //Transfer_matrix.mv(coarse->data(), fine->data());
        /*std::cout <<  "interpolate fein" <<  std::endl;
        for (int i=0; i< fine->data().size(); i++){
          std::cout <<  fine->data()[i] <<  std::endl;
        }*/
        std::cout <<  "+++ interpolate ende" <<  std::endl;        MPI_Barrier(MPI_COMM_WORLD);
	//std::exit(0);

        
      }
	std::cout <<  "+++ interpolate ende" <<  std::endl;
    }

    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().N();
      const size_t fine_ndofs = fine->get_data().N();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        // do a shortcut

	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();
	//if(my_rank==0) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine " << *r <<std::endl;} std::exit(0);
      } else {

	        std::cout <<  "+++ restrict " <<  std::endl;
   	/*std::cout <<  "restriktion fein" <<  std::endl;
        for (int i=0; i< fine->data().size(); i++){
          std::cout <<  fine->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction " <<  std::endl;*/
	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	if (my_rank==0) Dune::PDELab::Backend::native(fine->data())[0][0] = 0;
	if (my_rank==num_pro-1) Dune::PDELab::Backend::native(fine->data())[Dune::PDELab::Backend::native(fine->data()).N()-1][Dune::PDELab::Backend::native(fine->data()).N()-1] = 0;
        //std::cout <<  "ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg im restrict " <<  std::endl;
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank==1) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine " << *r <<std::endl;} 
	//inj->apply(fine->data(), coarse->data());
        inject_matrix.mv(Dune::PDELab::Backend::native(fine->data()), Dune::PDELab::Backend::native(coarse->data()));
	MPI_Barrier(MPI_COMM_WORLD); 
	if(my_rank==1) for(auto r =coarse->data().begin(); r !=coarse->data().end(); ++r){std::cout << "coarse " << *r <<std::endl;}// std::exit(0);
	//restrict_matrix.mtv(fine->data(), coarse->data());
        //Transfer_matrix2.mtv(fine->data(), coarse->data());
        //coarse->data() *= 0.5;
    /*std::cout <<  "restriction grob" <<  std::endl;
        for (int i=0; i< coarse->data().size(); i++){
          std::cout <<  coarse->data()[i] <<  std::endl;
        }
        std::cout <<  "restriction ende " <<  std::endl;*/


      }
    }
    
    template<class TransferTraits>
    void
    SpectralTransfer<
      TransferTraits>::restrict_u(const shared_ptr<typename TransferTraits::fine_encap_t> fine,
                                       shared_ptr<typename TransferTraits::coarse_encap_t> coarse)
    {
      ML_CVLOG(1, "TRANS", "restrict data");

      const size_t coarse_ndofs = coarse->get_data().N();
      const size_t fine_ndofs = fine->get_data().N();
      assert(coarse_ndofs > 0);
      assert(fine_ndofs >= coarse_ndofs);

      if (fine_ndofs == coarse_ndofs) {
        ML_CLOG(DEBUG, "TRANS", "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT");
        coarse->data() = fine->get_data();
	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	if(my_rank==0) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine " << *r <<std::endl;} //std::exit(0);
      } else {

	//rest->apply(fine->data(), coarse->data());
	//interpolate_matrix.mtv(fine->data(), coarse->data());
        //restrict_matrix.mv(Dune::PDELab::Backend::native(fine->data()), Dune::PDELab::Backend::native(coarse->data()));
	//if(my_rank==0) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine " << *r <<std::endl;} 	
        std::cout <<  "+++ restrict_u " <<  std::endl;

	int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );
	if (my_rank==0) Dune::PDELab::Backend::native(fine->data())[0][0] = 0;
	if (my_rank==num_pro-1) Dune::PDELab::Backend::native(fine->data())[Dune::PDELab::Backend::native(fine->data()).N()-1][Dune::PDELab::Backend::native(fine->data()).N()-1] = 0;

	if(my_rank==0) for(auto r =fine->data().begin(); r !=fine->data().end(); ++r){std::cout << "fine " << *r <<std::endl;} 



	for(int i =1; i<  Dune::PDELab::Backend::native(fine->data()).N()-1; i++){
			if(i==1){
						Dune::PDELab::Backend::native(fine->data())[i] = Dune::PDELab::Backend::native(coarse->data())[i*2-1] + 0.5 *Dune::PDELab::Backend::native(coarse->data())[i*2];}
			else if(i == Dune::PDELab::Backend::native(fine->data()).N()-1){} 
			else{Dune::PDELab::Backend::native(fine->data())[i] = 0.5*Dune::PDELab::Backend::native(coarse->data())[i*2-2] + Dune::PDELab::Backend::native(coarse->data())[i*2-1] + 0.5 *Dune::PDELab::Backend::native(coarse->data())[i*2];}
	}


	if(my_rank==0) for(auto r =coarse->data().begin(); r !=coarse->data().end(); ++r){std::cout << "coarse " << *r <<std::endl;} //std::exit(0);
	//if(my_rank==0) for(auto r =coarse->data().begin(); r !=coarse->data().end(); ++r){std::cout << "coarse " << *r <<std::endl;} //std::exit(0);

        


      }
    }
    
    
    
    
    
    
  }  // ::pfasst::contrib
}  // ::pfasst



