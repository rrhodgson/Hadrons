/*
 * DMixingUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MContraction_DMixingUtils_hpp_
#define Hadrons_MContraction_DMixingUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MContraction)

template <typename FImpl>
class DMixingUtils
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    /*     OpStruct
     * \   /      \    /
     *  \ /        \  /
     *   *          **
     *   *         /  \
     *  / \       /    \
     * /   \    
     *  One        Two
    */
    enum OpStruct {
        One=0,
        Two=1
    };
    static int toInt(OpStruct r);
    enum Parity {
        Pos=0,
        Neg=1
    };
    static std::string toString(Parity p);

    static const std::array<Gamma, 8> GHs;
    static const std::array<Gamma, 2> parityG;
    static void GH_cap(
        std::vector<PropagatorField> &out,
        const PropagatorField &prop,
        const Parity p);
    static void GH_cap(
        PropagatorField &out,
        const PropagatorField &prop,
        const Parity p,
        const OpStruct r);
};


template <typename FImpl>
int DMixingUtils<FImpl>::toInt(OpStruct r) {
    if (r == OpStruct::One)
        return 1;
    else if (r == OpStruct::Two)
        return 2;
    else
        HADRONS_ERROR(Argument, "DMixingUtils: Invalid OpStruct value");
}

template <typename FImpl>
const std::array<Gamma, 8> DMixingUtils<FImpl>::GHs = {
    Gamma(Gamma::Algebra::GammaX),
    Gamma(Gamma::Algebra::GammaY),
    Gamma(Gamma::Algebra::GammaZ),
    Gamma(Gamma::Algebra::GammaT),
    Gamma(Gamma::Algebra::GammaXGamma5),
    Gamma(Gamma::Algebra::GammaYGamma5),
    Gamma(Gamma::Algebra::GammaZGamma5),
    Gamma(Gamma::Algebra::GammaTGamma5)
};

template <typename FImpl>
const std::array<Gamma, 2> DMixingUtils<FImpl>::parityG = {
    Gamma(Gamma::Algebra::Identity),
    Gamma(Gamma::Algebra::Gamma5)
};

template <typename FImpl>
std::string DMixingUtils<FImpl>::toString(Parity p) {
    if (p == Parity::Pos)
        return "+";
    else if (p == Parity::Neg)
        return "-";
    else
        HADRONS_ERROR(Argument, "DMixingUtils: Invalid Parity value");
}


template <typename FImpl>
void DMixingUtils<FImpl>::GH_cap(
    std::vector<PropagatorField> &out,
    const PropagatorField &prop,
    const Parity p)
{
    assert(out.size() == 2);

    SitePropagator spId(1.0);
    out[OpStruct::One] = Zero();
    out[OpStruct::Two] = Zero();

    for (const auto &GH : GHs)
    {
        out[OpStruct::One] += GH * prop * parityG[p] * GH;
        out[OpStruct::Two] += spId * GH * trace(prop * parityG[p] * GH);
    }
};

template <typename FImpl>
void DMixingUtils<FImpl>::GH_cap(
    PropagatorField &out,
    const PropagatorField &prop,
    const Parity p,
    const OpStruct r)
{
    SitePropagator spId(1.0);
    out = Zero();

    switch (r)
    {
        case OpStruct::One:
            for (const auto &GH : GHs)
                out += GH * prop * parityG[p] * GH;
            break;

        case OpStruct::Two:
            for (const auto &GH : GHs)
                out += spId * GH * trace(prop * parityG[p] * GH);
            break;

        default:
            HADRONS_ERROR(Argument, "DMixingUtils: Invalid OpStruct value");
    }
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
