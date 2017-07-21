#include "FE_sweeper.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      template<
        class SweeperTrait,
        class BaseFunction,
        typename Enabled = void
      >
      class test
        : public Heat_FE<SweeperTrait, BaseFunction, Enabled>{
            
            
        public:
            explicit test<SweeperTrait, BaseFunction, Enabled>(std::shared_ptr<Dune::Functions::PQkNodalBasis<GridType::LevelGridView,SweeperTrait::BASE_ORDER>> basis, size_t nlevel, std::shared_ptr<GridType> grid)
        :   Heat_FE<SweeperTrait, BaseFunction, Enabled>(basis, nlevel, grid){}
            
        
          test(const test<SweeperTrait, BaseFunction, Enabled>& other) = default;
          
          test(test<SweeperTrait, BaseFunction, Enabled>&& other) = default;
          
          virtual ~test() = default;
            
            
            
            
        };   
    }
  }
}
