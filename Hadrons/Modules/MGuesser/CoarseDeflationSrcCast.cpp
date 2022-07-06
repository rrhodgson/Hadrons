#include <Hadrons/Modules/MGuesser/CoarseDeflationSrcCast.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TCoarseDeflationSrcCast<FIMPL,ARG(CoarseFermionEigenPack<FIMPLF, HADRONS_DEFAULT_LANCZOS_NBASIS>)>;
template class Grid::Hadrons::MGuesser::TCoarseDeflationSrcCast<FIMPL,ARG(CoarseFermionEigenPack<FIMPLF, 30>)>;
template class Grid::Hadrons::MGuesser::TCoarseDeflationSrcCast<FIMPL,ARG(CoarseFermionEigenPack<FIMPLF, 150>)>;
template class Grid::Hadrons::MGuesser::TCoarseDeflationSrcCast<FIMPL,ARG(CoarseFermionEigenPack<FIMPLF, 250>)>;
template class Grid::Hadrons::MGuesser::TCoarseDeflationSrcCast<FIMPL,ARG(CoarseFermionEigenPack<FIMPLF, 400>)>;
