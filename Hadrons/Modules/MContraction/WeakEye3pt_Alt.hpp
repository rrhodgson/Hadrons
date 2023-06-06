/*
 * WeakEye3pt_Alt.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: ferben <ferben@debian.felix.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MContraction_WeakEye3pt_Alt_hpp_
#define Hadrons_MContraction_WeakEye3pt_Alt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 * Weak Hamiltonian meson 3-pt diagrams, non-eye topologies.
 * 
 * Schematic:     
 *            qbl           qbr            |            qbl             qbr
 *          /--<--¬       /--<--¬          |          /--<--¬         /--<--¬       
 *         /       \     /       \         |         /       \       /       \      
 *        /         \   /         \        |        /         \     /         \     
 *       /           \ /           \       |       /           \   /           \    
 *  gIn *             * G           * gOut |  gIn *           G * * G           * gOut
 *      \             * G           |      |       \           /   \           /
 *       \           / \           /       |        \         /     \         /    
 *        \         /   \         /        |         \       /       \       /  
 *         \       /     \       /         |          \-->--/         \-->--/      
 *          \-->--/       \-->--/          |            ql              qr 
 *            ql            qr             |
 *               one trace                 |              two traces
 *
 * one trace : tr(ql*adj(gIn)*g5*adj(qbl)*g5*G*qbr*gOut*g5*adj(qr)*g5*G)
 * two traces: tr(ql*adj(gIn)*g5*adj(qbl)*g5*G)*tr(qbr*gOut*g5*adj(qr)*g5*G)
 * 
 */

BEGIN_MODULE_NAMESPACE(MContraction)

class WeakEye3pt_AltPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WeakEye3pt_AltPar,
                                    std::string,    q_loop,
                                    std::string,    qs_left,
                                    std::string,    qd_right,
                                    std::string,    qu,
                                    unsigned int,   ti,
                                    std::string,    xH,
                                    Gamma::Algebra, gammaIn,
                                    Gamma::Algebra, gammaOut,
                                    std::string,    output);
};

template <typename FImpl>
class TWeakEye3pt_Alt: public Module<WeakEye3pt_AltPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        Gamma::Algebra, in,
                                        Gamma::Algebra, out,
                                        Gamma::Algebra, GammaH1,
                                        Gamma::Algebra, GammaH2,
                                        unsigned int,   trace,
                                        unsigned int,   ti);
    };
    typedef Correlator<Metadata> Result;
public:
    // constructor
    TWeakEye3pt_Alt(const std::string name);
    // destructor
    virtual ~TWeakEye3pt_Alt(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WeakEye3pt_Alt, TWeakEye3pt_Alt<FIMPL>, MContraction);

/******************************************************************************
 *                 TWeakEye3pt_Alt implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWeakEye3pt_Alt<FImpl>::TWeakEye3pt_Alt(const std::string name)
: Module<WeakEye3pt_AltPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWeakEye3pt_Alt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q_loop, par().qs_left, par().qd_right, par().qu};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWeakEye3pt_Alt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TWeakEye3pt_Alt<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt_Alt<FImpl>::setup(void)
{
    envTmpLat(ComplexField, "corr");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWeakEye3pt_Alt<FImpl>::execute(void)
{
    LOG(Message) << "Computing mesonic weak 3pt contractions, non-eye topologies" << std::endl;
    LOG(Message) << "gIn : " << par().gammaIn  << std::endl;
    LOG(Message) << "gOut: " << par().gammaOut << std::endl;
    LOG(Message) << "q_loop  : " << par().q_loop   << std::endl;
    LOG(Message) << "qs_left : " << par().qs_left  << std::endl;
    LOG(Message) << "qd_right: " << par().qd_right << std::endl;
    LOG(Message) << "qu      : " << par().qu       << std::endl;

    std::vector<Result> result;
    Result              r;
    auto                &q   = envGet(PropagatorField, par().q_loop);
    auto                &qsl = envGet(PropagatorField, par().qs_left);
    auto                &qdr = envGet(PropagatorField, par().qd_right);
    auto                &qu  = envGet(PropagatorField, par().qu);


    SitePropagator      q_pt;
    SlicedPropagator    qsl_slice;
    SlicedPropagator    qdr_slice;
    SlicedPropagator    qu_slice;

    std::vector<int> xH = strToVec<int>(par().xH);

    peekSite(q_pt,q,xH);
    sliceSum(qsl,qsl_slice,Tp);
    sliceSum(qdr,qdr_slice,Tp);
    sliceSum(qu,qu_slice,Tp);

    Gamma               gIn(par().gammaIn), gOut(par().gammaOut);
    Gamma               g5(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, corr);
    r.info.in  = par().gammaIn;
    r.info.out = par().gammaOut;
    
    const std::array<std::pair<const Gamma::Algebra,const Gamma::Algebra>, 16> gH = {{
      { Gamma::Algebra::GammaX       , Gamma::Algebra::GammaX       },
      { Gamma::Algebra::GammaY       , Gamma::Algebra::GammaY       },
      { Gamma::Algebra::GammaZ       , Gamma::Algebra::GammaZ       },
      { Gamma::Algebra::GammaT       , Gamma::Algebra::GammaT       },
      { Gamma::Algebra::GammaX       , Gamma::Algebra::GammaXGamma5 },
      { Gamma::Algebra::GammaY       , Gamma::Algebra::GammaYGamma5 },
      { Gamma::Algebra::GammaZ       , Gamma::Algebra::GammaZGamma5 },
      { Gamma::Algebra::GammaT       , Gamma::Algebra::GammaTGamma5 },
      { Gamma::Algebra::GammaXGamma5 , Gamma::Algebra::GammaX       },
      { Gamma::Algebra::GammaYGamma5 , Gamma::Algebra::GammaY       },
      { Gamma::Algebra::GammaZGamma5 , Gamma::Algebra::GammaZ       },
      { Gamma::Algebra::GammaTGamma5 , Gamma::Algebra::GammaT       },
      { Gamma::Algebra::GammaXGamma5 , Gamma::Algebra::GammaXGamma5 },
      { Gamma::Algebra::GammaYGamma5 , Gamma::Algebra::GammaYGamma5 },
      { Gamma::Algebra::GammaZGamma5 , Gamma::Algebra::GammaZGamma5 },
      { Gamma::Algebra::GammaTGamma5 , Gamma::Algebra::GammaTGamma5 }
    }};

    int nt = qsl_slice.size();
    int ti = par().ti;

    for (auto& GH : gH)
    {
        const Gamma& GH1 = Gamma(GH.first);
        const Gamma& GH2 = Gamma(GH.second);

        r.info.GammaH1 = GH1.g;
        r.info.GammaH2 = GH2.g;

        {
            // one trace
            r.corr.clear();
            r.info.ti = ti;
            for (unsigned int tf=0; tf<nt; ++tf)
            {
                r.corr.push_back( TensorRemove(trace( GH1*q_pt*GH2*g5*adj(qdr_slice[tf])*g5*gOut*qu_slice[tf]*gIn*qsl_slice[ti])) );
            }
            r.info.trace = 1;
            result.push_back(r);
        }

        {
            // two trace
            r.corr.clear();
            r.info.ti = ti;
            for (unsigned int tf=0; tf<nt; ++tf)
            {
                r.corr.push_back( TensorRemove(trace( GH1*q_pt)*trace(GH2*g5*adj(qdr_slice[tf])*g5*gOut*qu_slice[tf]*gIn*qsl_slice[ti])) );
            }
            r.info.trace = 2;
            result.push_back(r);
        }

    }
    saveResult(par().output, "weakEye3pt_Alt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_WeakEye3pt_Alt_hpp_
