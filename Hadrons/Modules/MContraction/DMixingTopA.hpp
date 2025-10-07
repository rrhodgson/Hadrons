/*
 * WeakEye3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@cern.ch>
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
#ifndef Hadrons_MContraction_DMixingTopA_hpp_
#define Hadrons_MContraction_DMixingTopA_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopA                                        *
 *                (Fig. 4 (A) in arxiv:2504.16189)  
 *                   qCL        qInt1        qUR
 *                 /--<--\   /----<----\   /--<--\
 *                /       \ /           \ /       \
 *               /       ┌───┐         ┌───┐       \
 *           g5 *        | r |         | r'|        * g5
 *               \       └───┘         └───┘       /
 *                \       / \           / \       /
 *                 \-->--/   \---->----/   \-->--/
 *                   qUL        qInt1        qCR
 *
 * four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 * rr'=11: tr()
 * rr'=12: tr()
 * rr'=21: tr()
 * rr'=22: tr()
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopAPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopAPar,
                                    std::string,    qULeft,
                                    std::string,    qCLeft,
                                    std::string,    qURight,
                                    std::string,    qCRight,
                                    std::string,    qInt,
                                    std::string,    points,
                                    std::string,    output);
};

template <typename FImpl>
class TDMixingTopA: public Module<DMixingTopAPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, rr,
                                        std::string, parity);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;
public:
    // constructor
    TDMixingTopA(const std::string name);
    // destructor
    virtual ~TDMixingTopA(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual SlicedPropagator contractA_half_ti(const PropagatorField &GcuG, const std::vector<Coordinate>& xs); 
    virtual std::vector<SlicedPropagator> contractA_half_tf(const PropagatorField &GcuG, const std::vector<PropagatorField> &ds_prop_pt); 
    virtual std::vector<std::vector<Complex>> contractA(const SlicedPropagator &A, const std::vector<SlicedPropagator> &B); 
};

MODULE_REGISTER_TMP(DMixingTopA, TDMixingTopA<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopA implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopA<FImpl>::TDMixingTopA(const std::string name)
: Module<DMixingTopAPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft, 
	                           par().qCLeft,
	                           par().qURight,
	                           par().qCRight,
	                           par().qInt,
	                           par().points};
    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

template <typename FImpl>
typename TDMixingTopA<FImpl>::SlicedPropagator TDMixingTopA<FImpl>::contractA_half_ti(const TDMixingTopA<FImpl>::PropagatorField &GcuG, const std::vector<Coordinate>& xs) 
{
    int Nt = GcuG.Grid()->_fdimensions[3];
    SlicedPropagator A(Nt);
    for (int t1=0; t1<Nt; t1++) {
        A[t1] = peekSite(GcuG,xs[t1]);
    }
    return A;
}

template <typename FImpl>
std::vector<typename TDMixingTopA<FImpl>::SlicedPropagator> TDMixingTopA<FImpl>::contractA_half_tf(const TDMixingTopA<FImpl>::PropagatorField &GcuG, const std::vector<typename TDMixingTopA<FImpl>::PropagatorField> &ds_prop_pt) 
{
    int Nt = GcuG.Grid()->_fdimensions[3];
    Gamma g5(Gamma::Algebra::Gamma5);
    
    std::vector<SlicedPropagator> B(Nt, SlicedPropagator(Nt));
    for (int t1=0; t1<Nt; t1++) {
        SlicedPropagator buf;
        const auto& ds = ds_prop_pt[t1];
        PropagatorField tmp = g5*adj(ds)*g5 * GcuG * ds;
        sliceSum(tmp, buf, Tp);
        for (int t2=0; t2<Nt; t2++)
	{
            B[t1][t2] = buf[t2];
	}
    }
    return B;
}

template <typename FImpl>
std::vector<std::vector<Complex>> TDMixingTopA<FImpl>::contractA(const typename TDMixingTopA<FImpl>::SlicedPropagator &A, const std::vector<typename TDMixingTopA<FImpl>::SlicedPropagator> &B) 
{
    int Nt = A.size();
    
    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt));
    for (int t1=0; t1<Nt; t1++)
    {
        for (int t2=0; t2<Nt; t2++)
	{
            corr[t1][t2] = TensorRemove(trace( A[t1] * B[t1][t2] ));
	}
    }
    return corr;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::setup(void)
{
    envTmpLat(PropagatorField, "qcul");
    envTmpLat(PropagatorField, "qcur");
    envTmp(std::vector<PropagatorField>, "GcuG_l", 1, 2, PropagatorField(env().getGrid()));
    envTmp(std::vector<PropagatorField>, "GcuG_r", 1, 2, PropagatorField(env().getGrid()));
    const int Nt{env().getDim(Tdir)};
    envTmp(std::vector<SlicedPropagator>, "half_l", 1, 2, SlicedPropagator(Nt));
    envTmp(std::vector<std::vector<SlicedPropagator>>, "half_r", 1, 2, std::vector<SlicedPropagator>(Nt, SlicedPropagator(Nt)));

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qInt    : " << par().qInt << std::endl;
    LOG(Message) << "points  : " << par().points << std::endl;
      
    std::vector<Result> result;
    Result res;
    
    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);
    
    auto &qul    = envGet(PropagatorField, par().qULeft);
    auto &qcl    = envGet(PropagatorField, par().qCLeft);
    auto &qur    = envGet(PropagatorField, par().qURight);
    auto &qcr    = envGet(PropagatorField, par().qCRight);
    auto &qi     = envGet(std::vector<PropagatorField>, par().qInt);
    auto &points = envGet(std::vector<Coordinate>, par().points);

    envGetTmp(PropagatorField, qcul);
    envGetTmp(PropagatorField, qcur);
    envGetTmp(std::vector<PropagatorField>, GcuG_l);
    envGetTmp(std::vector<PropagatorField>, GcuG_r);
    envGetTmp(std::vector<SlicedPropagator>, half_l);
    envGetTmp(std::vector<std::vector<SlicedPropagator>>, half_r);

    Gamma g5(Gamma::Algebra::Gamma5);   

    qcul = qcl * g5 * g5 * adj(qul) * g5;
    qcur = qcr * g5 * g5 * adj(qur) * g5;

    DMixingUtils<FImpl>::GH_VVAA_cap(qcul, GcuG_l);
    DMixingUtils<FImpl>::GH_VVAA_cap(qcur, GcuG_r);

    std::vector<std::vector<Complex>> tmpRes = std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt, 0.));

    // parity +, parity -
    std::vector<Gamma> parityG = {Gamma(Gamma::Algebra::Identity), Gamma(Gamma::Algebra::Gamma5)};

    for (int p = 0; p < 2; p++)
    {
        for(int r = 0; r < 2; r++)
	{
            half_l[r] = contractA_half_ti(GcuG_l[r] * parityG[p], points);
            half_r[r] = contractA_half_tf(GcuG_r[r] * parityG[p], qi); 
	}
        for(int r = 0; r < 2; r++)
	{
            for(int s = 0; s < 2; s++)
            {
                res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);
                res.info.parity = (p == 0) ? "+" : "-";

	        tmpRes = contractA(half_l[r], half_r[s]);

		res.corr.clear();
                res.corr = tmpRes;
                result.push_back(res);

            }
        }
    }
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopA_hpp_
