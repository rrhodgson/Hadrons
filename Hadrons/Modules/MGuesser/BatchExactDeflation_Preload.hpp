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
                                    unsigned int, epSize,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename FImpl, typename EPack>
class TBatchExactDeflation_Preload: public Module<BatchExactDeflation_PreloadPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   PackField;
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

MODULE_REGISTER_TMP(BatchExactDeflation_Preload, ARG(TBatchExactDeflation_Preload<FIMPL,BaseFermionEigenPack<FIMPL>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflation_PreloadF, ARG(TBatchExactDeflation_Preload<FIMPLF,BaseFermionEigenPack<FIMPLF>>), MGuesser);
MODULE_REGISTER_TMP(BatchExactDeflation_Preload_OuterF, ARG(TBatchExactDeflation_Preload<FIMPL,BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl, typename EPack>
class BatchExactDeflation_PreloadGuesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field PackField;
public:
    BatchExactDeflation_PreloadGuesser(const std::vector<PackField> & evec,const std::vector<RealD> & eval, 
                                        const unsigned int epSize,
                                        const unsigned int evBatchSize,
                                        const unsigned int sourceBatchSize,
                                        const unsigned int traj)
    : evec_(evec)
    , eval_(eval)
    , epSize_(epSize)
    , evBatchSize_(evBatchSize)
    , sourceBatchSize_(sourceBatchSize)
    , traj_(traj)
    {};

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        unsigned int nBatch = epSize_/evBatchSize_ + (((epSize_ % evBatchSize_) != 0) ? 1 : 0);
        unsigned int sourceSize = out.size();

        std::vector<Field> evec_cast(evBatchSize_, Field(in[0].Grid()) );
        std::vector<RealD> eval_cast(evBatchSize_, 0.);

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        LOG(Message) << "--- zero guesses" << std::endl;
        for (auto &v: out)
        {
            v = Zero();
        }

        double cast_t = 0.;
        double proj_t = 0.;

        for (unsigned int bv = 0; bv < epSize_; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epSize_ - bv, evBatchSize_);

            cast_t -= usecond();
            if (!std::is_same<Field,PackField>::value) {
                for (unsigned int i = 0; i < evBlockSize; ++i) {
                    precisionChange(evec_cast[i],evec_[bv+i]);
                    eval_cast[i] = eval_[bv+i];
                }
            }
            cast_t += usecond();

            proj_t -= usecond();
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                if (std::is_same<Field,PackField>::value) {
                    projAccumulate(in, out, evec_, eval_, bv, bv + evBlockSize, bs, bs + sourceBlockSize);
                } else {
                    projAccumulate(in, out, evec_cast, eval_cast, 0, evBlockSize, bs, bs + sourceBlockSize);
                }
            }
            proj_t += usecond();
        }

        std::cout << "Total precision change time " << cast_t << std::endl;
        std::cout << "Total projection time " << proj_t << std::endl;
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    void projAccumulateImpl(const std::vector<Field> &in, std::vector<Field> &out,
                        const std::vector<Field>& evec_cast,
                        const std::vector<RealD>& eval_cast,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf,
                        std::true_type)
    {
        GridBase *g       = in[0].Grid();
        double   lVol     = g->lSites();
        double   siteSize = sizeof(typename Field::scalar_object);
        double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
        double   nIt      = (ef-ei)*(sf - si);
        double   t        = 0.;

        t -= usecond();
        for (unsigned int i = ei; i < ef; ++i)
        for (unsigned int j = si; j < sf; ++j)
        {
            axpy(out[j], 
                 TensorRemove(innerProduct(evec_cast[i], in[j]))/eval_cast[i], 
                 evec_cast[i], out[j]);
        }
        t += usecond();
        // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
        LOG(Message) << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;
    };

    template<typename F1, typename F2>
    void projAccumulateImpl(const std::vector<F1> &in, std::vector<F1> &out,
                        const std::vector<F2>& evec_cast,
                        const std::vector<RealD>& eval_cast,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf,
                        std::false_type)
    {
        assert(0 && "Type mismatch");
    };


    template<typename F1, typename F2>
    void projAccumulate(const std::vector<F1> &in, std::vector<F1> &out,
                        const std::vector<F2>& evec_cast,
                        const std::vector<RealD>& eval_cast,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf) {
        projAccumulateImpl(in, out,
                        evec_cast,
                        eval_cast,
                        ei, ef,
                        si, sf,
                        std::integral_constant<bool, std::is_same<Field,PackField>::value>{});
    }
private:
    const std::vector<PackField> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int          evBatchSize_, sourceBatchSize_;
    unsigned int          epSize_;
    int                   traj_;
};


/******************************************************************************
 *                     TBatchExactDeflation_Preload implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TBatchExactDeflation_Preload<FImpl, EPack>::TBatchExactDeflation_Preload(const std::string name)
: Module<BatchExactDeflation_PreloadPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflation_Preload<FImpl, EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflation_Preload<FImpl, EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TBatchExactDeflation_Preload<FImpl, EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflation_Preload<FImpl, EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().epSize << " modes) and batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflation_PreloadGuesser<FImpl,EPack>),
                     getName(), env().getObjectLs(par().eigenPack),
                     epack.evec, epack.eval, par().epSize, par().evBatchSize,
                     par().sourceBatchSize, vm().getTrajectory());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflation_Preload<FImpl, EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflation_Preload_hpp_
