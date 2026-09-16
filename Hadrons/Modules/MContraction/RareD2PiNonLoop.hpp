/*
 * RareD2PiNonLoop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_RareD2PiNonLoop_hpp_
#define Hadrons_MContraction_RareD2PiNonLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         RareD2PiNonLoop                                        *
 *                (Fig. 4 (A) in arxiv:2504.16189)
 *                   qCL       qUR
 *                 /-->--\   /-->--\
 *                /       \ /       \
 *               /       ┌───┐       \
 *           g5 *        | r |        * g5
 *               \       └───┘       /
 *                \       / \       /
 *                 \--<--/   \--<--/
 *                   qUL       qCR
 *         tsrc           t1            t2          tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB1 * qInt2 * GA2 * qCR * g5 * qUR * GB2]
 * rr' = 12:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB1 * qInt2 * GA2]*[qCR * g5 * qUR * GB2]
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class RareD2PiNonLoopPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RareD2PiNonLoopPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, output);
};

template <typename FImpl>
class TRareD2PiNonLoop : public Module<RareD2PiNonLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        int        , r,
                                        std::string, parity);
    };
    typedef Correlator<Metadata, Complex> Result;

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

public:
    // constructor
    TRareD2PiNonLoop(const std::string name);
    // destructor
    virtual ~TRareD2PiNonLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual std::vector<Complex> contract(const PropagatorField &GcuGl, const PropagatorField &qcur);
};

MODULE_REGISTER_TMP(RareD2PiNonLoop, TRareD2PiNonLoop<FIMPL>, MContraction);

/******************************************************************************
 *                        TRareD2PiNonLoop implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRareD2PiNonLoop<FImpl>::TRareD2PiNonLoop(const std::string name)
    : Module<RareD2PiNonLoopPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRareD2PiNonLoop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft,
                                   par().qCLeft,
                                   par().qURight,
                                   par().qCRight};
    return in;
}

template <typename FImpl>
std::vector<std::string> TRareD2PiNonLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TRareD2PiNonLoop<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
std::vector<Complex> TRareD2PiNonLoop<FImpl>::contract(
    const PropagatorField &GcuG_l,
    const PropagatorField &qcur)
{
    int Nt = env().getDim(Tdir);
    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<Complex> corr(Nt);

    SlicedComplex B;
    B.reserve(Nt);

    startTimer("mult");
    ComplexField tmp = trace(GcuG_l * qcur);
    stopTimer("mult");

    startTimer("sliceSum");
    sliceSum(tmp, B, Tp);
    stopTimer("sliceSum");

    startTimer("trace");
    for (int t = 0; t < Nt; t++)
        corr[t] = TensorRemove(B[t]);
    stopTimer("trace");

    return corr;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareD2PiNonLoop<FImpl>::setup(void)
{
    envTmpLat(PropagatorField, "qcul");
    envTmpLat(PropagatorField, "qcur");
    envTmp(std::vector<PropagatorField>, "GcuG_l", 1, 2, PropagatorField(env().getGrid()));
    envTmpLat(PropagatorField, "half_l");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRareD2PiNonLoop<FImpl>::execute(void)
{
    LOG(Message) << "Computing D2Pi NonLoop diagram" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;

    std::vector<Result> result;

    const int Nt{env().getDim(Tdir)};

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);

    envGetTmp(PropagatorField, qcul);
    envGetTmp(PropagatorField, qcur);
    envGetTmp(std::vector<PropagatorField>, GcuG_l);
    envGetTmp(PropagatorField, half_l);

    Gamma g5(Gamma::Algebra::Gamma5);

    qcul = qcl * adj(qul) * g5; // qcl * g5 * (g5 * adj(qul) * g5)
    qcur = qcr * adj(qur) * g5; // qcr * g5 * (g5 * adj(qur) * g5)

    startTimer("GH_cap");
    DMixingUtils<FImpl>::GH_cap(GcuG_l, qcul, Parity::Pos);
    stopTimer("GH_cap");

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (const auto r : {OpStruct::One,OpStruct::Two})
        {
            half_l = parityG[p] * GcuG_l[r];

            Result res;
            res.info.parity = DMixingUtils<FImpl>::toString(p);
            res.info.r      = DMixingUtils<FImpl>::toInt(r);

            res.corr = contract(half_l, qcur);

            result.push_back(res);
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "RareD2PiNonLoop", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_RareD2PiNonLoop_hpp_
