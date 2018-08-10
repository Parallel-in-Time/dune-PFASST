#ifndef VVAC_FUNCTIONS_HH
#define VVAC_FUNCTIONS_HH

#include <random>
#include <functional>
#include <dune/common/fvector.hh>


// basic functions ************************************************************************************************


// Random Disc for the scalar case
template <int dim>
class ScalarRandomDiscFunction
{
public:

    using DomainType = Dune::FieldVector<double, dim>;
    using RangeType = Dune::FieldVector<double, 1>;

    ScalarRandomDiscFunction(int seed, int n, double minRadius, double maxRadius, double lower, double upper, double value=1.0) :
        value_(value)
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> randomComponent(0, 1);
        std::uniform_real_distribution<> randomRadius(minRadius, maxRadius);
        std::uniform_real_distribution<> randomCoordinate(lower, upper);

        discs_.reserve(n);
        component_.reserve(n);
        for (int i=0; i<n; ++i)
        {
            DomainType center;
            for (int j=0; j<dim; ++j)
                center[j] = randomCoordinate(gen);
            double radius = randomRadius(gen);

            discs_.push_back([=](DomainType x)->bool {
                x -= center;
                return x.two_norm()<=radius;
            });
            component_.push_back(randomComponent(gen));
        }
    }

    void evaluate(const DomainType& x, RangeType& y) const
    {
        y = 0;
        y[0] = value_;
        for(std::size_t i=0; i<discs_.size(); ++i)
        {
            if (discs_[i](x))
            {
                //y[0] = 0;
                //y[component_[i]] = value_;
                y[0]=-1.0 + 2.0*component_[i]; // -1 if component zero, 1.0 if component one
                return;
            }
        }
    }

protected:
    std::vector<std::function<bool(DomainType)> > discs_;
    std::vector<int> component_;
    double value_;
};

#endif
