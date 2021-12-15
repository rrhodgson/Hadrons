#ifndef Hadrons_MGuesser_BatchExactDeflationPreloadUpCast_hpp_
#define Hadrons_MGuesser_BatchExactDeflationPreloadUpCast_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                           Batch deflation module                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class BatchExactDeflationPreloadUpCastPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BatchExactDeflationPreloadUpCastPar,
                                    std::string, eigenPack,
                                    unsigned int, epSize,
                                    unsigned int, evBatchSize,
                                    unsigned int, sourceBatchSize);
};

template <typename FImpl, typename EPack>
class TBatchExactDeflationPreloadUpCast: public Module<BatchExactDeflationPreloadUpCastPar>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
public:
    // constructor
    TBatchExactDeflationPreloadUpCast(const std::string name);
    // destructor
    virtual ~TBatchExactDeflationPreloadUpCast(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BatchExactDeflationPreloadUpCast, ARG(TBatchExactDeflationPreloadUpCast<FIMPL,BaseFermionEigenPack<FIMPLF>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template <typename FImpl, typename EPack>
class BatchExactDeflationPreloadUpCastGuesser: public LinearFunction<typename FImpl::FermionField>
{
public:
    typedef typename FImpl::FermionField Field;
    typedef typename EPack::Field   EPackField;
private:
    static constexpr bool requireCast_ = !(std::is_same<Field,EPackField>::value);
    static constexpr std::size_t CastBufferSize_( std::size_t n ){ return requireCast_ ? n : 0; };

    template <typename O>
    typename std::enable_if<std::is_same<Field,O>::value, const std::vector<Field> &>::type
    CopyOrOriginal( const std::vector<Field> &Copy, const std::vector<O> &Original ) { return Original; }
    template <typename O>
    typename std::enable_if<!std::is_same<Field,O>::value, const std::vector<Field> &>::type
    CopyOrOriginal( const std::vector<Field> &Copy, const std::vector<O> &Original ) { return Copy; }
public:
    BatchExactDeflationPreloadUpCastGuesser(const std::vector<EPackField> & evec,const std::vector<RealD> & eval, 
                                        const unsigned int epSize,
                                        const unsigned int evBatchSize,
                                        const unsigned int sourceBatchSize)
    : evec_(evec)
    , eval_(eval)
    , epSize_(epSize)
    , evBatchSize_(evBatchSize)
    , sourceBatchSize_(sourceBatchSize)
    {};

    virtual void operator() (const Field &in, Field &out)
    {}

    virtual void operator() (const std::vector<Field> &in, std::vector<Field> &out)
    {
        unsigned int nBatch = epSize_/evBatchSize_ + (((epSize_ % evBatchSize_) != 0) ? 1 : 0);
        unsigned int sourceSize = out.size();

        // std::vector<Field> evecCast( CastBufferSize_(evBatchSize_) , Field(in[0].Grid()) );
        // std::vector<RealD> evalCast( CastBufferSize_(evBatchSize_) , 0. );

        LOG(Message) << "=== BATCH DEFLATION GUESSER START" << std::endl;
        // if (requireCast_) {
        //     LOG(Message) << "Eigenpack requires precision change" << std::endl;
        // }

        for (auto &v: out)
        {
            v = Zero();
        }

        double cast_t = 0.;
        double proj_t = 0.;


        std::vector<EPackField> inCast( sourceSize , EPackField(evec_[0].Grid()) );
        std::vector<EPackField> outCast( sourceSize , EPackField(evec_[0].Grid()) );

        cast_t -= usecond();
        for (unsigned int i = 0; i < sourceSize; ++i) {
            precisionChange(inCast[i],in[i]);
            outCast[i] = Zero();
        }
        cast_t += usecond();

        for (unsigned int bv = 0; bv < epSize_; bv += evBatchSize_)
        {
            unsigned int evBlockSize = std::min(epSize_ - bv, evBatchSize_);

            proj_t -= usecond();
            for (unsigned int bs = 0; bs < sourceSize; bs += sourceBatchSize_)
            {
                unsigned int sourceBlockSize = std::min(sourceSize - bs, sourceBatchSize_);

                projAccumulate(inCast, outCast, 
                            evec_, eval_,
                            bv, bv + evBlockSize,
                            bs, bs + sourceBlockSize);
                // projAccumulate(inCast, outCast, 
                //             evec_, eval_,
                //             bv, bv + evBlockSize,
                //             0, sourceBlockSize);

            }
            proj_t += usecond();
        }

        cast_t -= usecond();
        for (unsigned int i = 0; i < sourceSize; ++i) {
            precisionChange(out[i],outCast[i]);
        }
        cast_t += usecond();

        LOG(Message) << "Total precision change time " << cast_t/1.e6 << " s" << std::endl;
        LOG(Message) << "Total projection time " << proj_t/1.e6 << " s" <<  std::endl;
        
        LOG(Message) << "=== BATCH DEFLATION GUESSER END" << std::endl;
    }
private:
    static void projAccumulate(const std::vector<EPackField> &in, std::vector<EPackField> &out,
                        const std::vector<EPackField>& evec,
                        const std::vector<RealD>& eval,
                        const unsigned int ei, const unsigned int ef,
                        const unsigned int si, const unsigned int sf)
    {
        GridBase *g       = in[0].Grid();
        double   lVol     = g->lSites();
        double   siteSize = sizeof(typename EPackField::scalar_object);
        double   lSizeGB  = lVol*siteSize/1024./1024./1024.;
        double   nIt      = (ef-ei)*(sf - si);
        double   t        = 0.;

        t -= usecond();
        for (unsigned int i = ei; i < ef; ++i)
        for (unsigned int j = si; j < sf; ++j)
        {
            axpy(out[j], 
                 TensorRemove(innerProduct(evec[i], in[j]))/eval[i], 
                 evec[i], out[j]);
        }
        t += usecond();
        // performance (STREAM convention): innerProduct 2 reads + axpy 2 reads 1 write = 5 transfers
        LOG(Message) << "projAccumulate: " << t << " us | " << 5.*nIt*lSizeGB << " GB | " << 5.*nIt*lSizeGB/t*1.0e6 << " GB/s" << std::endl;
    };

private:
    const std::vector<EPackField> &  evec_;
    const std::vector<RealD> &  eval_;
    unsigned int epSize_, evBatchSize_, sourceBatchSize_;
};


/******************************************************************************
 *                     TBatchExactDeflationPreloadUpCast implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
TBatchExactDeflationPreloadUpCast<FImpl, EPack>::TBatchExactDeflationPreloadUpCast(const std::string name)
: Module<BatchExactDeflationPreloadUpCastPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationPreloadUpCast<FImpl, EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};
    
    return in;
}

template <typename FImpl, typename EPack>
std::vector<std::string> TBatchExactDeflationPreloadUpCast<FImpl, EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl, typename EPack>
DependencyMap TBatchExactDeflationPreloadUpCast<FImpl, EPack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationPreloadUpCast<FImpl, EPack>::setup(void)
{
    LOG(Message) << "Setting batch exact deflation guesser with preloaded eigenPack '" << par().eigenPack << "'"
                 << "' (" << par().epSize << " modes) and batch size " 
                 << par().evBatchSize << ", and source batch size " 
                 << par().sourceBatchSize << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);

    envCreateDerived(LinearFunction<Field>, ARG(BatchExactDeflationPreloadUpCastGuesser<FImpl,EPack>),
                     getName(), env().getObjectLs(par().eigenPack),
                     epack.evec, epack.eval, par().epSize, par().evBatchSize,
                     par().sourceBatchSize);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename EPack>
void TBatchExactDeflationPreloadUpCast<FImpl, EPack>::execute(void)
{
    
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_BatchExactDeflationPreloadUpCast_hpp_
