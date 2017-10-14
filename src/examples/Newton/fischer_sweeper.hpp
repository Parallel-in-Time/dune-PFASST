#include "FE_sweeper.hpp"


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
          
          shared_ptr<typename SweeperTrait::encap_t> evaluate_rhs_impl(const typename SweeperTrait::time_t& m,const shared_ptr<typename SweeperTrait::encap_t> u) {
	
            ML_CVLOG(4, this->get_logger_id(),  "evaluating IMPLICIT part at t=" << m);
                    std::cout << "evaluate start" << std::endl;

            	auto result = this->get_encap_factory().create();
		auto newton_rhs = this->get_encap_factory().create();
		auto f = this->get_encap_factory().create();

          	result->zero();


		auto u_old = this->get_encap_factory().create();
		for(int k=0; k< this->last_newton_state()[0][m]->data().size(); k++){
    		u_old->data()[k] = this->last_newton_state()[0][m]->data()[k];}

		Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
            	evaluate_f2(f, u_old);
            	evaluate_df2(df, u_old);
            	df.mv(u_old->data(), newton_rhs->data());
            	newton_rhs->data() -= f->data();
		df.mv(u->data(), result->data());
		result->data() *= -1;
		result->data() -= newton_rhs->data();



		auto neu = this->get_encap_factory().create();
		neu->zero();
		this->df_dune[0][m]->mv(u->data(), neu->data());
		neu->data() += this->coarse_rhs()[0][m]->data(); 
		neu->data() *=-1;


		for (int i=0; i<u->get_data().size(); ++i)
          	{
     	     		std::cout << "vorher " << neu->data()[i] <<std::endl;	 
          	}

          	/*for (int i=0; i<u->get_data().size(); ++i)
          	{
     	     		fneu2[i] = pow(u->get_data()[i], _n + 1);	 
          	}
          	this->M_dune.mv(fneu2, result->data());
          
          	this->M_dune.mmv(u->get_data(), result->data());

          	result->data() *= (_nu*_nu);

	  	this->A_dune.mmv(u->get_data(),result->data());	

	  	result->data()*=-1;
        	for (size_t i = 0; i < u->get_data().size(); i++) {
	  		std::cout << "f evaluate " << result->data()[i] << std::endl;
        	}

                std::cout << "evaluate ende fast" << std::endl;
         	auto dt = this->get_status()->get_dt();
                std::cout << "evaluate ende dt bestimmt" << std::endl;

		for (int i=0; i<u->get_data().size(); ++i)
          	{
     	     		std::cout << result->get_data()[i] <<std::endl;	 
          	}*/

		/*evaluate_f(f, u, dt, rhs);
                evaluate_df(df, u, dt);
                df.mv(u->data(), newton_rhs);
                newton_rhs -= f->data();*/

      		/*if (this->is_coarse){
			std::cout << "grob" << std::endl;
        		result->data() +=  this->_M_initial->get_data(); //   this->get_states().front()->get_data();
          
      		}else{
        		this->M_dune.umv(this->get_states().front()->get_data(), result->data());
      		}

		std::cout << "vor _q delta" << std::endl;
      		for (size_t n = 0; n < (int) m; ++n) { //<=
			std::cout << "vor _q delta" << std::endl;
        		result->scaled_add(dt * this->_q_delta_impl(m + 1, n), this->_impl_rhs[n]);
      		}*/

                std::cout << "evaluate ende" << std::endl;     


            	return neu; //result


        }
         
        void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                                    shared_ptr<typename SweeperTrait::encap_t> u,
                                                    const typename SweeperTrait::time_t& m,
                                                    const typename SweeperTrait::time_t& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_t> rhs){
	
         ML_CVLOG(4, this->get_logger_id(), "IMPLICIT spatial SOLVE at node number" << m << " with dt=" << dt);
        std::cout << "implcit solve start" << std::endl;

         auto residuum = this->get_encap_factory().create();
	 Dune::BlockVector<Dune::FieldVector<double,1> > newton_rhs, newton_rhs2 ;
         newton_rhs.resize(rhs->get_data().size());
         newton_rhs2.resize(rhs->get_data().size());
    
	 //u->data() = this->last_newton_state()[0][m]->data();

    	//for(int i=0; i< num_time_steps; i++){	

		for(int k=0; k< this->last_newton_state()[0][m]->data().size(); k++){
		//std::cout << " " << m << " " << u->data()[k] << "   " << this->last_newton_state()[0][m]->data()[k] << std::endl;
    		u->data()[k] = this->last_newton_state()[0][m]->data()[k];}

	//}


         //u->zero();
	 for (int i=0; i< 1 ;i++){
            Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(this->M_dune); ///////M
            evaluate_f(f, u, dt, rhs);
            evaluate_df(df, u, dt);
            df.mv(u->data(), newton_rhs);
            newton_rhs -= f->data();
            newton_rhs2 = newton_rhs;

		
	    auto nv = this->get_encap_factory().create();
	    nv->data() = this->coarse_rhs()[0][m]->data();
	    nv->data() *= -dt;
	    nv->data() += rhs->data();

	    Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > df_neu = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >(*this->df_dune[0][m]); ///////M
 	    df_neu*=dt;
            df_neu+=this->M_dune;
            for (int i=0; i<u->get_data().size(); ++i)
          	{
     	     		//std::cout << "newton " << newton_rhs[i] << " " << nv->data()[i] <<std::endl;	 
          	}
            Dune::MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(df_neu); //neu df
	  
            Dune::SeqILU0<MatrixType,VectorType,VectorType> preconditioner(df_neu,1.0); //neu df
	  
            Dune::CGSolver<VectorType> cg(linearOperator,
                              preconditioner,
                              1e-16, // desired residual reduction factor
                              5000,    // maximum number of iterations
                              0);    // verbosity of the solver
          
          
            Dune::InverseOperatorResult statistics ;
            cg.apply(u->data(), nv->data() , statistics ); //newton_rhs

            evaluate_f(f, u, dt, rhs);
          
            //std::cout << i << " residuumsnorm von f(u) " << f->norm0() << std::endl;  
            //if(f->norm0()<1e-10){   break;} //  std::exit(0); std::cout << "genauigkeit erreicht " << i << std::endl;
          
            //df.mv(u->data(), residuum->data());
            //residuum->data() -= newton_rhs2;
            //std::cout << "residuums norm " << residuum->norm0() << std::endl;

	}

	Dune::BlockVector<Dune::FieldVector<double,1> > M_u;
        M_u.resize(u->get_data().size());
	this->M_dune.mv(u->get_data(), M_u);

        for (size_t i = 0; i < u->get_data().size(); i++) {
          f->data()[i] = (M_u[i] - rhs->get_data()[i]) / (dt);
	  //std::cout << "f " << f->data()[i] << std::endl;
        }

        //f = this->evaluate_rhs_impl(dt * this->_q_delta_impl(m+1, m+1), u); 
        //evaluate_rhs_impl(0, u);
	//std::exit(0);
        this->_num_impl_solves++;
                std::cout << "implcit solve end" << std::endl;

      }
      
        //private:
            
      void evaluate_f(shared_ptr<typename SweeperTrait::encap_t> f,
            const shared_ptr<typename SweeperTrait::encap_t> u,
            const typename SweeperTrait::time_t& dt,
            const shared_ptr<typename SweeperTrait::encap_t> rhs){

          double _nu=this->_nu;
          double _n=this->_n;

          f->zero();
          Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > fneu(this->M_dune);
          Dune::BlockVector<Dune::FieldVector<double,1> > fneu2;
          fneu2.resize(u->get_data().size());

          for (int i=0; i<u->get_data().size(); ++i)
          {
     	     fneu2[i] = pow(u->get_data()[i], _n + 1);	 
          }
          this->M_dune.mv(fneu2, f->data());
          
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
		for(int j=0; j< df.M(); j++)
                	if (df.exists(i,j)) df[i][j]= (_nu*_nu)*(_n+1) * pow(u->get_data()[j], _n) * this->M_dune[i][j];	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;
            df*=dt;
            df+=this->M_dune;
      } 
        


      void evaluate_f2(shared_ptr<typename SweeperTrait::encap_t> f,
            const shared_ptr<typename SweeperTrait::encap_t> u){

          double _nu=this->_nu;
          double _n=this->_n;

          f->zero();
          Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > fneu(this->M_dune);
          Dune::BlockVector<Dune::FieldVector<double,1> > fneu2;
          fneu2.resize(u->get_data().size());

          for (int i=0; i<u->get_data().size(); ++i)
          {
     	     fneu2[i] = pow(u->get_data()[i], _n + 1);	 
          }
          this->M_dune.mv(fneu2, f->data());
          
          this->M_dune.mmv(u->get_data(), f->data());

          f->data() *= (_nu*_nu);


          this->A_dune.mmv(u->get_data(),f->data());
          
      }
						
      void evaluate_df2(Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > &df,
                                                 const shared_ptr<typename SweeperTrait::encap_t> u
 						){
            double _nu=this->_nu;
            double _n=this->_n;
            for (int i=0; i<df.N(); ++i)
            {
		for(int j=0; j< df.M(); j++)
                	if (df.exists(i,j)) df[i][j]= (_nu*_nu)*(_n+1) * pow(u->get_data()[j], _n) * this->M_dune[i][j];	
            }
            df.axpy((-_nu*_nu), this->M_dune);
            df-=this->A_dune;

      } 






  
    };   
    }
  }
}
