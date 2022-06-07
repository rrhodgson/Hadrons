#ifndef Hadrons_MGuesser_ExactDeflationSrcCast_hpp_
#define Hadrons_MGuesser_ExactDeflationSrcCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Exact deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class ExactDeflationSrcCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDeflationSrcCastPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename FImpl, typename EPack>
class TExactDeflationSrcCast: public Module<ExactDeflationSrcCastPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field        FieldF;
public:
    // constructor
    TExactDeflationSrcCast(const std::string name);
    // destructor
    virtual ~TExactDeflationSrcCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDeflationSrcCast, ARG(TExactDeflationSrcCast<FIMPL, BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                 TExactDeflationSrcCast implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TExactDeflationSrcCast<FImpl,EPack>::TExactDeflationSrcCast(const std::string name)
: Module<ExactDeflationSrcCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TExactDeflationSrcCast<FImpl,EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TExactDeflationSrcCast<FImpl,EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TExactDeflationSrcCast<FImpl,EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}


template<class Field, class FieldF>
class DeflatedSrcCastGuesser: public LinearFunction<Field> {
private:
  const std::vector<FieldF> &evec;
  const std::vector<RealD> &eval;
  const unsigned int       N;

public:
  using LinearFunction<Field>::operator();

  DeflatedSrcCastGuesser(const std::vector<FieldF> & _evec,const std::vector<RealD> & _eval)
  : DeflatedSrcCastGuesser(_evec, _eval, _evec.size())
  {}

  DeflatedSrcCastGuesser(const std::vector<FieldF> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
  : evec(_evec), eval(_eval), N(_N)
  {
    assert(evec.size()==eval.size());
    assert(N <= evec.size());
  } 

  virtual void operator()(const Field &src,Field &guess) {
    FieldF src_f(evec[0].Grid());
    precisionChange(src_f,src);
    FieldF guess_f(src_f.Grid());
    guess_f = Zero();
    
    for (int i=0;i<N;i++) {
      const FieldF& tmp = evec[i];
      axpy(guess_f,TensorRemove(innerProduct(tmp,src_f)) / eval[i],tmp,guess_f);
    }
    guess_f.Checkerboard() = src_f.Checkerboard();

    // guess = guess_f; // prec change
    precisionChange(guess,guess_f);
  }
};


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TExactDeflationSrcCast<FImpl,EPack>::setup(void)
{
    LOG(Message) << "Setting exact deflation src cast guesser with eigenpack '"
                 << par().eigenPack << "' (" 
                 << par().size << " modes)" << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ARG(DeflatedSrcCastGuesser<Field,FieldF>), getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TExactDeflationSrcCast<FImpl,EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_ExactDeflationSrcCast_hpp_
