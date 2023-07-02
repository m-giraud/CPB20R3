// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
 /*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief A helper class to fill the flux variable caches used in the flux constitutive laws
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Helper class to fill the flux variables caches within
 *        the interaction volume around a given sub-control volume face.
 */
template<class TypeTag>
class CCMpfaFluxVariablesCacheFiller
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using MpfaHelper = typename FVGridGeometry::MpfaHelper;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;

    using PrimaryInteractionVolume = GetPropType<TypeTag, Properties::PrimaryInteractionVolume>;
    using PrimaryDataHandle = typename ElementFluxVariablesCache::PrimaryIvDataHandle;
    using PrimaryLocalFaceData = typename PrimaryInteractionVolume::Traits::LocalFaceData;
    using SecondaryInteractionVolume = GetPropType<TypeTag, Properties::SecondaryInteractionVolume>;
    using SecondaryDataHandle = typename ElementFluxVariablesCache::SecondaryIvDataHandle;
    using SecondaryLocalFaceData = typename SecondaryInteractionVolume::Traits::LocalFaceData;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static constexpr bool doAdvection = ModelTraits::enableAdvection();
    static constexpr bool doDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool doHeatConduction = ModelTraits::enableEnergyBalance();

    static constexpr bool soldependentAdvection = getPropValue<TypeTag, Properties::SolutionDependentAdvection>();
    static constexpr bool soldependentDiffusion = getPropValue<TypeTag, Properties::SolutionDependentMolecularDiffusion>();
    static constexpr bool soldependentHeatConduction = getPropValue<TypeTag, Properties::SolutionDependentHeatConduction>();

public:
    //! This cache filler is always solution-dependent, as it updates the
    //! vectors of cell unknowns with which the transmissibilities have to be
    //! multiplied in order to obtain the fluxes.
    static constexpr bool isSolDependent = true;

    //! The constructor. Sets problem pointer.
    CCMpfaFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheStorage Class that holds the scvf flux vars caches
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param ivDataStorage Class that stores the interaction volumes & handles
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables (primary/secondary variables)
     * \param scvf The corresponding sub-control volume face
     * \param forceUpdateAll if true, forces all caches to be updated (even the solution-independent ones)
     */
    template<class FluxVarsCacheStorage, class IVDataStorage>
    void fill(FluxVarsCacheStorage& fluxVarsCacheStorage,
              FluxVariablesCache& scvfFluxVarsCache,
              IVDataStorage& ivDataStorage,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              bool forceUpdateAll = false)
    {
        // Set pointers
        elementPtr_ = &element;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;

        // prepare interaction volume and fill caches of all the scvfs connected to it
        const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
        if (fvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = ivDataStorage.secondaryInteractionVolumes.size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().secondaryIndexSet(scvf);
                ivDataStorage.secondaryInteractionVolumes.emplace_back();
                secondaryIv_ = &ivDataStorage.secondaryInteractionVolumes.back();
                secondaryIv_->bind(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                ivDataStorage.secondaryDataHandles.emplace_back();
                secondaryIvDataHandle_ = &ivDataStorage.secondaryDataHandles.back();

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheStorage, *secondaryIv_, *secondaryIvDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                secondaryIv_ = &ivDataStorage.secondaryInteractionVolumes[ivIndexInContainer];
                secondaryIvDataHandle_ = &ivDataStorage.secondaryDataHandles[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheStorage, *secondaryIv_, *secondaryIvDataHandle_, ivIndexInContainer);
            }
        }
        else
        {
            if (forceUpdateAll)
            {
                // the local index of the interaction volume to be created in its container
                const auto ivIndexInContainer = ivDataStorage.primaryInteractionVolumes.size();

                // prepare the locally cached boundary interaction volume
                const auto& indexSet = fvGridGeometry.gridInteractionVolumeIndexSets().primaryIndexSet(scvf);
                ivDataStorage.primaryInteractionVolumes.emplace_back();
                primaryIv_ = &ivDataStorage.primaryInteractionVolumes.back();
                primaryIv_->bind(indexSet, problem(), fvGeometry);

                // prepare the corresponding data handle
                ivDataStorage.primaryDataHandles.emplace_back();
                primaryIvDataHandle_ = &ivDataStorage.primaryDataHandles.back();

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheStorage, *primaryIv_, *primaryIvDataHandle_, ivIndexInContainer, true);
            }
            else
            {
                const auto ivIndexInContainer = scvfFluxVarsCache.ivIndexInContainer();
                primaryIv_ = &ivDataStorage.primaryInteractionVolumes[ivIndexInContainer];
                primaryIvDataHandle_ = &ivDataStorage.primaryDataHandles[ivIndexInContainer];

                // fill the caches for all the scvfs in the interaction volume
                fillCachesInInteractionVolume_(fluxVarsCacheStorage, *primaryIv_, *primaryIvDataHandle_, ivIndexInContainer);
            }
        }
    }

    //! returns the stored interaction volume pointer
    const PrimaryInteractionVolume& primaryInteractionVolume() const
    { return *primaryIv_; }

    //! returns the stored interaction volume pointer
    const SecondaryInteractionVolume& secondaryInteractionVolume() const
    { return *secondaryIv_; }

    //! returns the stored data handle pointer
    const PrimaryDataHandle& primaryIvDataHandle() const
    { return *primaryIvDataHandle_; }

    //! returns the stored data handle pointer
    const SecondaryDataHandle& secondaryIvDataHandle() const
    { return *secondaryIvDataHandle_; }

    //! returns the currently stored iv-local face data object
    const PrimaryLocalFaceData& primaryIvLocalFaceData() const
    { return *primaryLocalFaceData_; }

    //! returns the currently stored iv-local face data object
    const SecondaryLocalFaceData& secondaryIvLocalFaceData() const
    { return *secondaryLocalFaceData_; }

private:

    const Problem& problem() const { return *problemPtr_; }
    const Element& element() const { return *elementPtr_; }
    const FVElementGeometry& fvGeometry() const { return *fvGeometryPtr_; }
    const ElementVolumeVariables& elemVolVars() const { return *elemVolVarsPtr_; }

    //! Method to fill the flux var caches within an interaction volume
    template<class FluxVarsCacheStorage, class InteractionVolume, class DataHandle>
    void fillCachesInInteractionVolume_(FluxVarsCacheStorage& fluxVarsCacheStorage,
                                        InteractionVolume& iv,
                                        DataHandle& handle,
                                        unsigned int ivIndexInContainer,
                                        bool forceUpdateAll = false)
    {
        // determine if secondary interaction volumes are used here
        static constexpr bool isSecondary = MpfaHelper::considerSecondaryIVs()
                                            && std::is_same<InteractionVolume, SecondaryInteractionVolume>::value;

        // First we upate data which are not dependent on the physical processes.
        // We store pointers to the other flux var caches, so that we have to obtain
        // this data only once and can use it again in the sub-cache fillers.
        const auto numGlobalScvfs = iv.localFaceData().size();
        std::vector<const SubControlVolumeFace*> ivScvfs(numGlobalScvfs);
        std::vector<FluxVariablesCache*> ivFluxVarCaches(numGlobalScvfs);

        unsigned int i = 0;
        for (const auto& d : iv.localFaceData())
        {
            // obtain the scvf
            const auto& scvfJ = fvGeometry().scvf(d.gridScvfIndex());
            ivScvfs[i] = &scvfJ;
            ivFluxVarCaches[i] = &fluxVarsCacheStorage[scvfJ];
            ivFluxVarCaches[i]->setIvIndexInContainer(ivIndexInContainer);
            ivFluxVarCaches[i]->setUpdateStatus(true);
            ivFluxVarCaches[i]->setSecondaryIvUsage(isSecondary);
            ivFluxVarCaches[i]->setIvLocalFaceIndex(d.ivLocalScvfIndex());
            if (dim < dimWorld)
                if (d.isOutsideFace())
                    ivFluxVarCaches[i]->setIndexInOutsideFaces(d.scvfLocalOutsideScvfIndex());
            i++;
        }

        fillAdvection(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
        fillDiffusion(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
        fillHeatConduction(iv, handle, ivScvfs, ivFluxVarCaches, forceUpdateAll);
    }

    //! fills the advective quantities (enabled advection)
    template< class InteractionVolume,
              class DataHandle,
              bool enableAdvection = doAdvection,
              typename std::enable_if_t<enableAdvection, int> = 0 >
    void fillAdvection(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {
        using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
        using AdvectionFiller = typename AdvectionType::Cache::Filler;

        // fill data in the handle
        fillAdvectionHandle(iv, handle, forceUpdateAll);

        // fill advection caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
            else
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

            // fill this scvfs cache
            AdvectionFiller::fill(*ivFluxVarCaches[i],
                                  problem(),
                                  iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                  fvGeometry(),
                                  elemVolVars(),
                                  *ivScvfs[i],
                                  *this);
        }
    }

    //! do nothing if advection is not enabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableAdvection = doAdvection,
              typename std::enable_if_t<!enableAdvection, int> = 0 >
    void fillAdvection(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {}

    //! fills the diffusive quantities (diffusion enabled)
    template< class InteractionVolume,
              class DataHandle,
              bool enableDiffusion = doDiffusion,
              typename std::enable_if_t<enableDiffusion, int> = 0 >
    void fillDiffusion(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {
        using DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
        using DiffusionFiller = typename DiffusionType::Cache::Filler;

        static constexpr int numPhases = ModelTraits::numFluidPhases();
        static constexpr int numComponents = ModelTraits::numFluidComponents();

        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            handle.diffusionHandle().setPhaseIndex(phaseIdx);
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
                if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                    continue;

                // fill data in the handle
                handle.diffusionHandle().setComponentIndex(compIdx);
                fillDiffusionHandle(iv, handle, forceUpdateAll, phaseIdx, compIdx);

                // fill diffusion caches
                for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
                {
                    // set pointer to current local face data object
                    // ifs are evaluated at compile time and are optimized away
                    if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
                    {
                        // we cannot make a disctinction, thus we set both pointers
                        primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                        secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
                    }
                    else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
                        primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                    else
                        secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

                    // fill this scvfs cache
                    DiffusionFiller::fill(*ivFluxVarCaches[i],
                                          phaseIdx,
                                          compIdx,
                                          problem(),
                                          iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                          fvGeometry(),
                                          elemVolVars(),
                                          *ivScvfs[i],
                                          *this);
                }
            }
        }
    }

    //! do nothing if diffusion is not enabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableDiffusion = doDiffusion,
              typename std::enable_if_t<!enableDiffusion, int> = 0 >
    void fillDiffusion(InteractionVolume& iv,
                       DataHandle& handle,
                       const std::vector<const SubControlVolumeFace*>& ivScvfs,
                       const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                       bool forceUpdateAll = false)
    {}

    //! fills the quantities related to heat conduction (heat conduction enabled)
    template< class InteractionVolume,
              class DataHandle,
              bool enableHeatConduction = doHeatConduction,
              typename std::enable_if_t<enableHeatConduction, int> = 0 >
    void fillHeatConduction(InteractionVolume& iv,
                            DataHandle& handle,
                            const std::vector<const SubControlVolumeFace*>& ivScvfs,
                            const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                            bool forceUpdateAll = false)
    {
        using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;
        using HeatConductionFiller = typename HeatConductionType::Cache::Filler;

        // prepare data in handle
        fillHeatConductionHandle(iv, handle, forceUpdateAll);

        // fill heat conduction caches
        for (unsigned int i = 0; i < iv.localFaceData().size(); ++i)
        {
            // set pointer to current local face data object
            // ifs are evaluated at compile time and are optimized away
            if (std::is_same<PrimaryInteractionVolume, SecondaryInteractionVolume>::value)
            {
                // we cannot make a disctinction, thus we set both pointers
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);
            }
            else if (std::is_same<InteractionVolume, PrimaryInteractionVolume>::value)
                primaryLocalFaceData_ = &(iv.localFaceData()[i]);
            else
                secondaryLocalFaceData_ = &(iv.localFaceData()[i]);

            // fill this scvfs cache
            HeatConductionFiller::fill(*ivFluxVarCaches[i],
                                       problem(),
                                       iv.element(iv.localFaceData()[i].ivLocalInsideScvIndex()),
                                       fvGeometry(),
                                       elemVolVars(),
                                       *ivScvfs[i],
                                       *this);
        }
    }

    //! do nothing if heat conduction is disabled
    template< class InteractionVolume,
              class DataHandle,
              bool enableHeatConduction = doHeatConduction,
              typename std::enable_if_t<!enableHeatConduction, int> = 0 >
    void fillHeatConduction(InteractionVolume& iv,
                            DataHandle& handle,
                            const std::vector<const SubControlVolumeFace*>& ivScvfs,
                            const std::vector<FluxVariablesCache*>& ivFluxVarCaches,
                            bool forceUpdateAll = false)
    {}

    //! prepares the quantities necessary for advective fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>,
              typename std::enable_if_t<AdvectionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillAdvectionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;

        // get instance of the interaction volume-local assembler
        using Traits = typename InteractionVolume::Traits;
        using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // Assemble T only if permeability is sol-dependent or if update is forced
        if (forceUpdateAll || soldependentAdvection)
            localAssembler.assembleMatrices(handle.advectionHandle(), iv, LambdaFactory::getAdvectionLambda());

        // assemble pressure vectors
        for (unsigned int pIdx = 0; pIdx < ModelTraits::numFluidPhases(); ++pIdx)
        {
            // set context in handle
            handle.advectionHandle().setPhaseIndex(pIdx);

            // maybe (re-)assemble gravity contribution vector
            auto getRho = [pIdx] (const auto& volVars) { return volVars.density(pIdx); };
            static const bool enableGravity = getParamFromGroup<bool>(problem().paramGroup(), "Problem.EnableGravity");
            if (enableGravity)
                localAssembler.assembleGravity(handle.advectionHandle(), iv, getRho);

            // reassemble pressure vector
            auto getPressure = [pIdx] (const auto& volVars) { return volVars.pressure(pIdx); };
            localAssembler.assembleU(handle.advectionHandle(), iv, getPressure);
        }
    }

    //! prepares the quantities necessary for diffusive fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>,
              typename std::enable_if_t<DiffusionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillDiffusionHandle(InteractionVolume& iv,
                             DataHandle& handle,
                             bool forceUpdateAll,
                             int phaseIdx, int compIdx)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;

        // get instance of the interaction volume-local assembler
        using Traits = typename InteractionVolume::Traits;
        using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // maybe (re-)assemble matrices
        if (forceUpdateAll || soldependentDiffusion)
            localAssembler.assembleMatrices(handle.diffusionHandle(),
                                            iv,
                                            LambdaFactory::getDiffusionLambda(phaseIdx, compIdx));

        // assemble vector of mole fractions
        auto getMoleFraction = [phaseIdx, compIdx] (const auto& volVars) { return volVars.moleFraction(phaseIdx, compIdx); };
        localAssembler.assembleU(handle.diffusionHandle(), iv, getMoleFraction);
    }

    //! prepares the quantities necessary for conductive fluxes in the handle
    template< class InteractionVolume,
              class DataHandle,
              class HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>,
              typename std::enable_if_t<HeatConductionType::discMethod == DiscretizationMethod::ccmpfa, int> = 0 >
    void fillHeatConductionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll)
    {
        using LambdaFactory = TensorLambdaFactory<DiscretizationMethod::ccmpfa>;
        using ThermCondModel = GetPropType<TypeTag, Properties::ThermalConductivityModel>;

        // get instance of the interaction volume-local assembler
        using Traits = typename InteractionVolume::Traits;
        using IvLocalAssembler = typename Traits::template LocalAssembler<Problem, FVElementGeometry, ElementVolumeVariables>;
        IvLocalAssembler localAssembler(problem(), fvGeometry(), elemVolVars());

        // maybe (re-)assemble matrices
        if (forceUpdateAll || soldependentHeatConduction)
            localAssembler.assembleMatrices(handle.heatConductionHandle(),
                                            iv,
                                            LambdaFactory::template getHeatConductionLambda<ThermCondModel>());

        // assemble vector of temperatures
        auto getMoleFraction = [] (const auto& volVars) { return volVars.temperature(); };
        localAssembler.assembleU(handle.heatConductionHandle(), iv, getMoleFraction);
    }

    //! fill handle only when advection uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>,
              typename std::enable_if_t<AdvectionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillAdvectionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll) {}

    //! fill handle only when diffusion uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class DiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>,
              typename std::enable_if_t<DiffusionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillDiffusionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll, int phaseIdx, int compIdx) {}

    //! fill handle only when heat conduction uses mpfa
    template< class InteractionVolume,
              class DataHandle,
              class HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>,
              typename std::enable_if_t<HeatConductionType::discMethod != DiscretizationMethod::ccmpfa, int> = 0 >
    void fillHeatConductionHandle(InteractionVolume& iv, DataHandle& handle, bool forceUpdateAll) {}

    const Problem* problemPtr_;
    const Element* elementPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    // We store pointers to an inner and a boundary interaction volume.
    // These are updated during the filling of the caches and the
    // physics-related caches have access to them
    PrimaryInteractionVolume* primaryIv_;
    SecondaryInteractionVolume* secondaryIv_;

    // pointer to the current interaction volume data handle
    PrimaryDataHandle* primaryIvDataHandle_;
    SecondaryDataHandle* secondaryIvDataHandle_;

    // We do an interaction volume-wise filling of the caches
    // While filling, we store a pointer to the current localScvf
    // face data object of the IV so that the individual caches
    // can access it and don't have to retrieve it again
    const PrimaryLocalFaceData* primaryLocalFaceData_;
    const SecondaryLocalFaceData* secondaryLocalFaceData_;
};

} // end namespace Dumux

#endif
