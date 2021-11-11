#ifndef Hadrons_MGuesser_Prop_Guesser_hpp_
#define Hadrons_MGuesser_Prop_Guesser_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MIO/LoadEigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class Prop_GuesserPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Prop_GuesserPar,
                                    std::string, prop);
};

template <typename FImpl>
class TProp_Guesser: public Module<Prop_GuesserPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
public:
    // constructor
    TProp_Guesser(const std::string name);
    // destructor
    virtual ~TProp_Guesser(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(PropGuesser, ARG(TProp_Guesser<FIMPL>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl>
class Prop_Guesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    Prop_Guesser(GridBase * grid, GridBase * gridRb)
    {
        grid_ = grid;
        gridRb_ = gridRb;
        fields.resize(12, gridRb_);

        // for (int i=0; i<12; i++) {
        //     fields[i] = Zero();
        // }

    };


    void setProp(PropagatorField& prop) {
        int j=0;
        FermionField tmp(grid_);

        for (int s=0; s<4; s++) {
        for (int c=0; c<3; c++) {
            std::cout << "Spin " << s << " , colour " << c << std::endl;
            // PropToFerm<FImpl>(fields[j], prop, s, c);

            PropToFerm<FImpl>(tmp, prop, s, c);
            // std::cout << s << " , " << c << std::endl << tmp << std::endl;
            pickCheckerboard(Odd ,fields[j],tmp);
            // std::cout << s << " , " << c << std::endl << fields[j] << std::endl;

            j++;
        }}
    }

    virtual void operator() (const FermionField &in, FermionField &out)
    {}

    virtual void operator() (const std::vector<FermionField> &in, std::vector<FermionField> &out)
    {
        LOG(Message) << "=== BATCH PROP GUESSER START" << std::endl;
        
        for (int i=0; i<out.size(); i++) {
            // fields[i].Checkerboard() = out[i].Checkerboard();
            out[i] = fields[i];
        }

        LOG(Message) << "=== BATCH PROP GUESSER END" << std::endl;
    }
private:
    GridBase * grid_;
    GridBase * gridRb_;
    std::vector<FermionField> fields;
};


/******************************************************************************
 *                     TProp_Guesser implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TProp_Guesser<FImpl>::TProp_Guesser(const std::string name)
: Module<Prop_GuesserPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TProp_Guesser<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    in.push_back(par().prop);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TProp_Guesser<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl>
DependencyMap TProp_Guesser<FImpl>::getObjectDependencies(void)
{
    DependencyMap dep;

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TProp_Guesser<FImpl>::setup(void)
{
    envCreateDerived(LinearFunction<FermionField>, Prop_Guesser<FImpl>,
                     getName(), 1, envGetGrid(FermionField, env().getObjectLs(par().prop) ), envGetRbGrid(FermionField, env().getObjectLs(par().prop) ));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TProp_Guesser<FImpl>::execute(void)
{

    // LOG(Message) << "Setting batch exact deflation guesser with eigenpack "
    //              << "located at '" << par().eigenPack.filestem
    //              << "' (" << par().eigenPack.size << " modes), EV batch size " 
    //              << par().evBatchSize << ", and source batch size " 
    //              << par().sourceBatchSize << std::endl;
    
    auto &prop      = envGet(PropagatorField, par().prop);

    auto &guesser   = envGetDerived(LinearFunction<FermionField>, Prop_Guesser<FImpl>, getName());

    guesser.setProp(prop);


    

    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_Prop_Guesser_hpp_