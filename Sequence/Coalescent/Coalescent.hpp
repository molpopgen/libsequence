#ifndef __SEQUENCE_COALESCENT_COALECENT_HPP_
#define __SEQUENCE_COALESCENT_COALECENT_HPP__

/*!
  \defgroup coalescent Classes and functions related to simulating data under coalescent models
  \ingroup popgen
*/
/*! \file Coalescent.hpp
  @brief A lazy header to include the headers needed to start writing simulations.
  Includes:
   <Sequence/Coalescent/SimTypes.hpp>
   <Sequence/Coalescent/Coalesce.hpp>
   <Sequence/Coalescent/Recombination.hpp>
   <Sequence/Coalescent/Mutation.hpp>
   <Sequence/Coalescent/Initialize.hpp>
   <Sequence/Coalescent/DemographicModels.hpp>
   <Sequence/Coalescent/FragmentsRescaling.hpp>
*/
/*! \example freerec.cc
  Coalescent simulation with free recombination
*/
/*! \example ms--.cc
  Coalescent simulation
*/
/*! \example msbeta.cc
  Coalescent simulation with non-uniform genetic maps
*/
/*! \example bottleneck.cc
  Example of using the  Sequence::bottleneck template function
*/
/*! \example fragments.cc
  Example of simulating partially linked fragments in neutral models.
*/

#include <Sequence/Coalescent/SimTypes.hpp>
#include <Sequence/Coalescent/Coalesce.hpp>
#include <Sequence/Coalescent/Recombination.hpp>
#include <Sequence/Coalescent/Mutation.hpp>
#include <Sequence/Coalescent/Initialize.hpp>
#include <Sequence/Coalescent/DemographicModels.hpp>
#include <Sequence/Coalescent/FragmentsRescaling.hpp>
#include <Sequence/Coalescent/Trajectories.hpp>
#endif
