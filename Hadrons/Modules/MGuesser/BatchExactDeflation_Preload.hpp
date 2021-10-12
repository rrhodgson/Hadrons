#ifndef Hadrons_MGuesser_BatchExactDeflation_Preload_hpp_
#define Hadrons_MGuesser_BatchExactDeflation_Preload_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflation_PreloadPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflation_PreloadPar,
                                    std::string, eigenPack,
                                    unsigned int, size,
                                    unsigned int, batchSize);
};

template <typename EPack>
class TBatchExactDeflation_Preload: public Module<BatchExactDeflation_PreloadPar>
{
public:
    typedef typename EPack::Field   Field;
public:
    // constructor
    TBatchExactDeflation_Preload(const std::string name);
    // destructor
    virtual ~TBatchExactDeflation_Preload(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflation_Preload, ARG(TBatchExactDeflation_Preload<BaseFermionEigenPack<FIMPL>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflation_PreloadF, ARG(TBatchExactDeflation_Preload<BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename EPack>
class BatchExactDeflation_PreloadGuesser: public LinearFunction<typename EPack::Field>
{
public:
    typedef typename EPack::Field Field;
public:
    BatchExactDeflation_PreloadGuesser(const std::vector<Field> & evec,const std::vector<RealD> & eval, 
                                        const unsigned int size,
                                        const unsigned int batchSize,
                                        const unsigned int traj)
    : evec_(evec)
    , eval_(eval)
    , size_(size)
    , batchSize_(batchSize)
    , traj_(traj)
    {

    };

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        unsigned int nBatch = size_/batchSize_ + (((size_ % batchSize_) != 0) ? 1 : 0);

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        LOG(Message) << "--- zero guesses" << std::endl;
        for (auto &v: out)
        {
            v = Zero();
        }
        for (unsigned int b = 0; b < size_; b += batchSize_)
        {
            unsigned int bsize = std::min(size_ - b, batchSize_);

            LOG(Message) << "--- batch " << b/batchSize_ << std::endl;
            LOG(Message) << "project" << std::endl;
            projAccumulate(in, out, b, bsize);
        }
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    void projAccumulate(const std::vector<Field> &in, std::vector<Field> &out,
                        const unsigned int startIdx, const unsigned int batchSize)
    {
        for (unsigned int i = 0; i < batchSize; ++i)
        for (unsigned int j = 0; j < in.size(); ++j)
        {
            unsigned int idx = startIdx + i;
            axpy(out[j], 
                 TensorRemove(innerProduct(evec_[idx], in[j]))/eval_[idx], 
                 evec_[idx], out[j]);
        }
    };
private:
    const std::vector<Field> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int          batchSize_;
    unsigned int          size_;
    int                   traj_;
};


/******************************************************************************
 *                     TBatchExactDeflation_Preload implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TBatchExactDeflation_Preload<EPack>::TBatchExactDeflation_Preload(const std::string name)
: Module<BatchExactDeflation_PreloadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TBatchExactDeflation_Preload<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename EPack>
std::vector<std::string> TBatchExactDeflation_Preload<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename EPack>
DependencyMap TBatchExactDeflation_Preload<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TBatchExactDeflation_Preload<EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().size << " modes) and batch size " 
                 << par().batchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflation_PreloadGuesser<EPack>),
                     getName(), env().getObjectLs(par().eigenPack), epack.evec, epack.eval, par().size, par().batchSize, vm().getTrajectory());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TBatchExactDeflation_Preload<EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflation_Preload_hpp_
