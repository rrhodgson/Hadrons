/*
 * RareD2PiLoop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2025
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matteo Di Carlo <matteo.dicarlo@cern.ch>
 * Author: Felix Erben <felix.erben@cern.ch>
 * Author: Raoul Hodgson <raoul.hodgson@desy.de>
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
#ifndef Hadrons_MContraction_RareD2PiLoop_hpp_
#define Hadrons_MContraction_RareD2PiLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)
 *             qCL   ┌───┐   qUR
 *         /---->----| r |---->----\
 *        /          └───┘          \
 *       /          /     \          \
 *   g5 *           \     /           * g5
 *       \           \->-/ qLoop    /
 *        \                         /
 *         \-----------<-----------/
 *                    qUL
 *  tsrc               t1              tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qCR * g5 * qUR * GB1 * qLoop * GA1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]
 * rr' = 12:
 *  [qCR * g5 * qUR * GB1 * qLoop * GA1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]
 *
 ******************************************************************************/

BEGIN_MODULE_NAMESPACE(MContraction)

class RareD2PiLoopPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RareD2PiLoopPar,
                                    std::string, qCLeft,
                                    std::string, qULeft,
                                    std::string, qURight,
                                    std::string, qLoop,
                                    unsigned int, tOut,
                                    std::string, output);
};

template <typename FImpl>
class TRareD2PiLoop : public Module<RareD2PiLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        int        , r,
                                        std::string, parity,
                                        int        , eta);
    };
    typedef Correlator<Metadata, Complex> Result;

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

public:
    // constructor
    TRareD2PiLoop(const std::string name);
    // destructor
    virtual ~TRareD2PiLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual std::vector<Complex> contract(
        const PropagatorField &prop_c,
        const SitePropagator  &prop_ul_adj,
        const PropagatorField &prop_ur_adj,
        const PropagatorField &loop);
};

MODULE_REGISTER_TMP(RareD2PiLoop, TRareD2PiLoop<FIMPL>, MContraction);

/******************************************************************************
 *                        TRareD2PiLoop implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRareD2PiLoop<FImpl>::TRareD2PiLoop(const std::string name)
    : Module<RareD2PiLoopPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRareD2PiLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qCLeft,
                                   par().qULeft,
                                   par().qURight};
    if ( !par().qLoop.empty() )
        in.push_back(par().qLoop);

    return in;
}

template <typename FImpl>
std::vector<std::string> TRareD2PiLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TRareD2PiLoop<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
std::vector<Complex> TRareD2PiLoop<FImpl>::contract(
    const PropagatorField &prop_c,
    const SitePropagator  &prop_ul_adj,
    const PropagatorField &prop_ur_adj,
    const PropagatorField &loop)
{
    Gamma g5(Gamma::Algebra::Gamma5);

    startTimer("trace");
    LatticeComplex tmp = trace(prop_c * g5 * prop_ul_adj * g5 * prop_ur_adj * loop);
    stopTimer("trace");

    startTimer("sliceSum");
    SlicedComplex ret;
    sliceSum(tmp, ret, Tp);
    stopTimer("sliceSum");

    const int Nt{env().getDim(Tdir)};
    std::vector<Complex> out(Nt);
    for (unsigned int t = 0; t < Nt; t++)
    {
        out[t] = TensorRemove(ret[t]);
    }
    return out;
};

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareD2PiLoop<FImpl>::setup(void)
{
    GridCartesian *grid = envGetGrid(FermionField);
    envTmp(std::vector<PropagatorField>, "GdsG", 1, 2, PropagatorField(env().getGrid()));

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareD2PiLoop<FImpl>::execute(void)
{
    LOG(Message) << "Computing RareD2Pi loop diagram" << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    if (par().qLoop.empty())
        LOG(Message) << "Empty qLoop : Scalar/Pseudoscalar bilinear" << std::endl;
    else
        LOG(Message) << "qLoop  : " << par().qLoop << std::endl;

    std::vector<Result> result;

    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);

    std::vector<PropagatorField *> ql_default(1,nullptr);

    const auto &qcl = envGet(PropagatorField, par().qCLeft);
    const auto &qul = envGet(SlicedPropagator, par().qULeft);
    const auto &qur = envGet(PropagatorField, par().qURight);
    auto &ql = (par().qLoop.empty())  ? ql_default
                                        : envGet(std::vector<PropagatorField *>, par().qLoop);

    PropagatorField Unit(grid); Unit = 1;

    Gamma g5(Gamma::Algebra::Gamma5);
    const SitePropagator  qul_adj = g5 * adj(qul[par().tOut]) * g5;
    const PropagatorField qur_adj = g5 * adj(qur) * g5;

    int Neta = ql.size();

    envGetTmp(std::vector<PropagatorField>, GdsG);

    std::vector<OpStruct> r_vals;
    if (par().qLoop.empty()) r_vals = {OpStruct::One};
    else                     r_vals = {OpStruct::One,OpStruct::Two};

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (int i = 0; i < Neta; i++)
        {
            startTimer("GH_cap");
            if (par().qLoop.empty()) {
                GdsG[OpStruct::One] = Unit * parityG[p];
                GdsG[OpStruct::Two] = Unit * parityG[p];
            } else {
                DMixingUtils<FImpl>::GH_cap(GdsG, *ql[i], p);
            }
            stopTimer("GH_cap");

            for (const auto r : r_vals)
            {
                Result res;
                res.info.parity = DMixingUtils<FImpl>::toString(p);
                res.info.r      = DMixingUtils<FImpl>::toInt(r);
                res.info.eta    = i;

                res.corr = contract(qcl, qul_adj, qur_adj, GdsG[r]);
                result.push_back(res);
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "RareD2PiLoop", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_RareD2PiLoop_hpp_
