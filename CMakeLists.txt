add_library(fabm_models_ogs OBJECT
            ogs_model_library.F90
            bfm_pelagic_base.F90
            Phyto.F90
            PelBac.F90
<<<<<<< HEAD
            PelChem.F90
=======
            MicroZoo.F90
>>>>>>> bdca1d0e79e5955febc548b09930acbf71c16f7d
            shared.F90
            light.F90
           )

add_dependencies(fabm_models_ogs fabm_base)
