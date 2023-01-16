#ifndef Hadrons_MUtilities_AddFields_hpp_
#define Hadrons_MUtilities_AddFields_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Pack fields as a vector of pointers                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class AddFieldsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AddFieldsPar,
                                    std::vector<std::string>, fields);
};

template <typename Field>
class TAddFields: public Module<AddFieldsPar>
{
public:
    // constructor
    TAddFields(const std::string name);
    // destructor
    virtual ~TAddFields(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(AddPropagatorFields, TAddFields<FIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                 TAddFields implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TAddFields<Field>::TAddFields(const std::string name)
: Module<AddFieldsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TAddFields<Field>::getInput(void)
{
    std::vector<std::string> in = par().fields;
    
    return in;
}

template <typename Field>
std::vector<std::string> TAddFields<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename Field>
DependencyMap TAddFields<Field>::getObjectDependencies(void)
{
    DependencyMap dep;

    for (auto &n: par().fields)
    {
        dep.insert({n, getName()});
    }

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TAddFields<Field>::setup(void)
{
    int Ls = 0;
    for (int i=0; i<par().fields.size(); i++) {
        int Ls_tmp = env().getObjectLs(par().fields[i]);
        if (i==0 || Ls == Ls_tmp) {
            Ls = Ls_tmp;
        } else {
            LOG(Message) << "i = " << i << "  :  " << Ls << " != " << Ls_tmp << std::endl;
            assert(0);
        }
    }
    envCreateLat(Field, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TAddFields<Field>::execute(void)
{
    LOG(Message) << "Adding " << par().fields.size() << " fields" << std::endl;
    LOG(Message) << "fields: " << par().fields << std::endl;
    auto &field = envGet(Field, getName());
    field = Zero();

    for (unsigned int i = 0; i < par().fields.size(); ++i)
    {
        auto& tmp = envGet(Field, par().fields[i]);
        field += tmp;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_AddFields_hpp_
