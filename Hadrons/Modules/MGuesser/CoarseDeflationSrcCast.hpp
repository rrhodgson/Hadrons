#ifndef Hadrons_MGuesser_CoarseDeflationSrcCast_hpp_
#define Hadrons_MGuesser_CoarseDeflationSrcCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Exact deflation guesser                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class CoarseDeflationSrcCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CoarseDeflationSrcCastPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename FImpl, typename EPack>
class TCoarseDeflationSrcCast: public Module<CoarseDeflationSrcCastPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field        FieldF;
    typedef typename EPack::CoarseField  CoarseFieldF;
public:
    // constructor
    TCoarseDeflationSrcCast(const std::string name);
    // destructor
    virtual ~TCoarseDeflationSrcCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

typedef CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS> CEPack;
typedef CoarseFermionEigenPack<FIMPLF,30> CEPack30;
typedef CoarseFermionEigenPack<FIMPLF,150> CEPack150;
typedef CoarseFermionEigenPack<FIMPLF,400> CEPack400;

MODULE_REGISTER_TMP(CoarseDeflationSrcCast, ARG(TCoarseDeflationSrcCast<FIMPL, CEPack>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflationSrcCast30, ARG(TCoarseDeflationSrcCast<FIMPL, CEPack30>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflationSrcCast150, ARG(TCoarseDeflationSrcCast<FIMPL, CEPack150>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflationSrcCast400, ARG(TCoarseDeflationSrcCast<FIMPL, CEPack400>), MGuesser);

/******************************************************************************
 *                 TCoarseDeflationSrcCast implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TCoarseDeflationSrcCast<FImpl,EPack>::TCoarseDeflationSrcCast(const std::string name)
: Module<CoarseDeflationSrcCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TCoarseDeflationSrcCast<FImpl,EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TCoarseDeflationSrcCast<FImpl,EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TCoarseDeflationSrcCast<FImpl,EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}


// template<class Field, class FieldF>
// class DeflatedSrcCastGuesser: public LinearFunction<Field> {
// private:
//   const std::vector<FieldF> &evec;
//   const std::vector<RealD> &eval;
//   const unsigned int       N;

// public:
//   using LinearFunction<Field>::operator();

//   DeflatedSrcCastGuesser(const std::vector<FieldF> & _evec,const std::vector<RealD> & _eval)
//   : DeflatedSrcCastGuesser(_evec, _eval, _evec.size())
//   {}

//   DeflatedSrcCastGuesser(const std::vector<FieldF> & _evec, const std::vector<RealD> & _eval, const unsigned int _N)
//   : evec(_evec), eval(_eval), N(_N)
//   {
//     assert(evec.size()==eval.size());
//     assert(N <= evec.size());
//   } 

//   virtual void operator()(const Field &src,Field &guess) {
//     FieldF src_f(evec[0].Grid());
//     precisionChange(src_f,src);
//     FieldF guess_f(src_f.Grid());
//     guess_f = Zero();
    
//     for (int i=0;i<N;i++) {
//       const FieldF& tmp = evec[i];
//       axpy(guess_f,TensorRemove(innerProduct(tmp,src_f)) / eval[i],tmp,guess_f);
//     }
//     guess_f.Checkerboard() = src_f.Checkerboard();

//     // guess = guess_f; // prec change
//     precisionChange(guess,guess_f);
//   }
// };

template<class FineField, class FineFieldF, class CoarseFieldF>
class LocalCoherenceDeflatedSrcCastGuesser: public LinearFunction<FineField> {
private:
  const std::vector<FineFieldF>   &subspace;
  const std::vector<CoarseFieldF> &evec_coarse;
  const std::vector<RealD>       &eval_coarse;
public:
  
  using LinearFunction<FineField>::operator();
  LocalCoherenceDeflatedSrcCastGuesser(const std::vector<FineFieldF>   &_subspace,
                const std::vector<CoarseFieldF> &_evec_coarse,
                const std::vector<RealD>       &_eval_coarse)
    : subspace(_subspace), 
      evec_coarse(_evec_coarse), 
      eval_coarse(_eval_coarse)  
  {
  }
  
  void operator()(const FineField &src,FineField &guess) { 
    std::cout << "Coarse deflating (src cast)" << std::endl;
    int N = (int)evec_coarse.size();

    FineFieldF src_f(subspace[0].Grid());
    precisionChange(src_f,src);
    FineFieldF guess_f(src_f.Grid());
    guess_f = Zero();

    CoarseFieldF src_coarse(evec_coarse[0].Grid());
    CoarseFieldF guess_coarse(evec_coarse[0].Grid());
    guess_coarse = Zero();

    blockProject(src_coarse,src_f,subspace);    
    for (int i=0;i<N;i++) {
      const CoarseFieldF & tmp = evec_coarse[i];
      axpy(guess_coarse,TensorRemove(innerProduct(tmp,src_coarse)) / eval_coarse[i],tmp,guess_coarse);
    }
    blockPromote(guess_coarse,guess_f,subspace);
    guess.Checkerboard() = src.Checkerboard();

    precisionChange(guess,guess_f);
  };
};


// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TCoarseDeflationSrcCast<FImpl,EPack>::setup(void)
{
    LOG(Message) << "Setting coarse deflation src cast guesser with eigenpack '"
                 << par().eigenPack << "' (" 
                 << par().size << " modes)" << std::endl;
    
    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ARG(LocalCoherenceDeflatedSrcCastGuesser<Field,FieldF,CoarseFieldF>), getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.evecCoarse, epack.evalCoarse);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TCoarseDeflationSrcCast<FImpl,EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_CoarseDeflationSrcCast_hpp_
