/*
 * DMixingTopD.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_DMixingTopD_hpp_
#define Hadrons_MContraction_DMixingTopD_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)
 *             qCL   ┌───┐                 qUR
 *         /---->----| r |----------------->--------------\
 *        /          └───┘                                 \
 *       /          /     \          qLoop2 /-<-\           \
 *   g5 *           \     /                /     \           * g5
 *       \           \->-/ qLoop1          \     /          /
 *        \                                 ┌───┐          /
 *         \--------------<-----------------| r'|----<----/
 *                       qUL                └───┘  qCR
 *  tsrc               t1                    t2              tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]
 * rr' = 12:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]
 * rr' = 21:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]*[qLoop1 * GA1]
 * rr' = 22:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]*[qLoop1 * GA1]
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopDPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopDPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, qLoop1,
                                    std::string, qLoop2,
                                    std::string, output);
};

template <typename FImpl>
class TDMixingTopD : public Module<DMixingTopDPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, rr,
                                        std::string, parity,
                                        std::string, eta_max);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

public:
    // constructor
    TDMixingTopD(const std::string name);
    // destructor
    virtual ~TDMixingTopD(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual SlicedPropagator contract_D_half(
        const PropagatorField &prop_c,
        const PropagatorField &prop_u_adj,
        const PropagatorField &loop);
    virtual std::vector<std::vector<Complex>> contract_D(
        const SlicedPropagator &half_if,
        const SlicedPropagator &half_fi);

    int half_idx(OpStruct r, int i, int Neta) {
        return i + Neta * r;
    }
};

MODULE_REGISTER_TMP(DMixingTopD, TDMixingTopD<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopD implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopD<FImpl>::TDMixingTopD(const std::string name)
    : Module<DMixingTopDPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft,
                                   par().qCLeft,
                                   par().qURight,
                                   par().qCRight,
                                   par().qLoop1,
                                   par().qLoop2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
typename TDMixingTopD<FImpl>::SlicedPropagator TDMixingTopD<FImpl>::contract_D_half(
    const PropagatorField &prop_c,
    const PropagatorField &prop_u_adj,
    const PropagatorField &loop)
{
    PropagatorField tmp = prop_u_adj * loop * prop_c;
    SlicedPropagator out;
    sliceSum(tmp, out, Tp);
    return out;
};

template <typename FImpl>
std::vector<std::vector<Complex>> TDMixingTopD<FImpl>::contract_D(
    const SlicedPropagator &half_l,
    const SlicedPropagator &half_r)
{
    Gamma g5(Gamma::Algebra::Gamma5);

    int Nt = env().getDim(Tdir);
    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt));

    for (int t1 = 0; t1 < Nt; t1++)
    {
        for (int t2 = 0; t2 < Nt; t2++)
        {
            corr[t1][t2] = TensorRemove(trace(half_l[t1] * g5 * half_r[t2] * g5));
        }
    }

    return corr;
};

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::setup(void)
{

    GridCartesian *grid = envGetGrid(FermionField);

    auto &ql1 = envGet(std::vector<PropagatorField *>, par().qLoop1);
    int Neta = ql1.size();
    const int Nt = env().getDim(Tdir);

    envTmp(std::vector<SlicedPropagator>, "half_l", 1, 2 * Neta, SlicedPropagator(Nt));
    envTmp(std::vector<SlicedPropagator>, "half_r", 1, 2 * Neta, SlicedPropagator(Nt));
    envTmp(std::vector<PropagatorField>, "GdsG", 1, 2, PropagatorField(env().getGrid()));

    if (par().qLoop1 != par().qLoop2)
    {
        HADRONS_ERROR(Argument, "Current implementation for identical loops only");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qLoop1  : " << par().qLoop1 << std::endl;
    LOG(Message) << "qLoop2  : " << par().qLoop2 << std::endl;

    std::vector<Result> result;

    const int Nt{env().getDim(Tdir)};

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &ql1 = envGet(std::vector<PropagatorField *>, par().qLoop1);
    auto &ql2 = envGet(std::vector<PropagatorField *>, par().qLoop2);

    Gamma g5(Gamma::Algebra::Gamma5);
    const PropagatorField qur_adj = g5 * adj(qur) * g5;
    const PropagatorField qul_adj = g5 * adj(qul) * g5;

    // this is assuming both loops are identical
    int Neta = ql1.size();
    bool same_loop_noise = true;

    envGetTmp(std::vector<SlicedPropagator>, half_l);
    envGetTmp(std::vector<SlicedPropagator>, half_r);
    envGetTmp(std::vector<PropagatorField>, GdsG);

    std::vector<std::vector<Complex>> tmpRes(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<std::vector<Complex>> tmpSum(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<std::vector<Complex>> diagSum(Nt, std::vector<Complex>(Nt, 0.));

    SlicedPropagator Lsum(Nt), Rsum(Nt);

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (int i = 0; i < Neta; i++)
        {
            // here one has to add ql2 if one wants them to be allowed to be different
            startTimer("GH_cap");
            DMixingUtils<FImpl>::GH_cap(GdsG, *ql1[i], p);
            stopTimer("GH_cap");
            for (const auto r : {OpStruct::One,OpStruct::Two})
            {
                startTimer("contract_D_half");
                half_l[half_idx(r,i,Neta)] = contract_D_half(qcl, qur_adj, GdsG[r]);
                half_r[half_idx(r,i,Neta)] = contract_D_half(qcr, qul_adj, GdsG[r]);
                stopTimer("contract_D_half");
            }
        }

        for (const auto r : {OpStruct::One,OpStruct::Two})
        {
            for (const auto s : {OpStruct::One,OpStruct::Two})
            {
                for (int t = 0; t < Nt; t++)
                {
                    Lsum[t] = Zero();
                    Rsum[t] = Zero();
                    if (same_loop_noise)
                        std::fill(diagSum[t].begin(), diagSum[t].end(), Complex(0.0));
                }

                for (int i = 0; i < Neta; i++) // imax = i + 1
                {
                    const auto &Li = half_l[half_idx(r,i,Neta)];
                    const auto &Ri = half_r[half_idx(s,i,Neta)];

                    // accumulate sums of Li, Ri
                    for (int t = 0; t < Nt; t++)
                    {
                        Lsum[t] += Li[t];
                        Rsum[t] += Ri[t];
                    }
                    startTimer("contract_D");
                    tmpSum = contract_D(Lsum, Rsum);
                    stopTimer("contract_D");

                    if (same_loop_noise)
                    {
                        // accumulate sum of diagonal terms & remove from total
                        startTimer("contract_D");
                        const auto diag = contract_D(Li, Ri);
                        stopTimer("contract_D");
                        for (int t1 = 0; t1 < Nt; t1++)
                            for (int t2 = 0; t2 < Nt; t2++)
                            {
                                diagSum[t1][t2] += diag[t1][t2];
                                tmpSum[t1][t2] -= diagSum[t1][t2];
                            }
                    }

                    const double norm = (!same_loop_noise)
                                            ? double((i + 1) * (i + 1))
                                            : (i > 0 ? double((i + 1) * i) : 1.0);

                    for (int t1 = 0; t1 < Nt; t1++)
                        for (int t2 = 0; t2 < Nt; t2++)
                            tmpRes[t1][t2] = tmpSum[t1][t2] / norm;

                    Result res;
                    res.info.parity = (p == Parity::Pos) ? "+" : "-";
                    res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);
                    res.info.eta_max = std::to_string(i + 1);
                    res.corr = tmpRes;
                    result.push_back(res);
                }
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopD", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopD_hpp_
