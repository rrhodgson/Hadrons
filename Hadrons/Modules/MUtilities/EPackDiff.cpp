/*
 * EPackDiff.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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
#include <Hadrons/Modules/MUtilities/EPackDiff.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MUtilities;

template class Grid::Hadrons::MUtilities::TEPackDiff<FermionEigenPack<FIMPL>, FIMPL, GIMPL>;
template class Grid::Hadrons::MUtilities::TEPackDiff<FermionEigenPack<FIMPLF>, FIMPLF, GIMPLF>;

// MODULE_REGISTER_TMP(EPackDiff, ARG(TEPackDiff<FermionEigenPack<FIMPL>, GIMPL>), MUtilities);
// MODULE_REGISTER_TMP(EPackDiffF, ARG(TEPackDiff<FermionEigenPack<FIMPLF>, GIMPLF>), MUtilities);