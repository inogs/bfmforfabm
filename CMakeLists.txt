add_library(fabm_models_ogs OBJECT
            ogs_model_library.F90
            bfm_pelagic_base.F90
            Phyto.F90
            PelBac.F90
            shared.F90
            light.F90
           )

add_dependencies(fabm_models_ogs fabm_base)
