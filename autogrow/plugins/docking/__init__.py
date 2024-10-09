from abc import abstractmethod
import glob
import os
import random
from typing import Any, Dict, List, Optional, Tuple, Union, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.docking.docking_class.parent_pdbqt_converter import ParentPDBQTConverter
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PostDockedCompoundInfo, PreDockedCompoundInfo, ScoreType
import autogrow.docking.scoring.execute_scoring_mol as Scoring
import autogrow.docking.ranking.ranking_mol as Ranking


class DockingBase(PluginBase):
    def run(self, **kwargs) -> Any:
        """Run the plugin with provided arguments."""
        lig_pdbqt_filename: str = kwargs["lig_pdbqt_filename"]
        file_conversion_class_object: ParentPDBQTConverter = kwargs[
            "file_conversion_class_object"
        ]

        return self.run_docking(
            lig_pdbqt_filename=lig_pdbqt_filename,
            file_conversion_class_object=file_conversion_class_object,
        )

    @abstractmethod
    def run_docking(
        self,
        lig_pdbqt_filename: str,
        file_conversion_class_object: ParentPDBQTConverter,
    ) -> Optional[str]:
        """
        run_docking is needs to be implemented in each class.

        Inputs:
        :param str pdbqt_filename: a string for docking process raise exception
            if missing

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        # raise NotImplementedError("run_dock() not implemented")
        pass


class DockingPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> Any:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: bool: True if the molecule passes the filter, False if it fails
        """

        dockings = self.get_selected_plugins_from_params()

        if dockings is None or len(dockings) == 0:
            raise Exception(
                f"You must specify a docking program! Choose from {str(self.plugins.keys())}"
            )
        if len(dockings) > 1:
            raise Exception(
                f"Only one docking program can be selected at a time! You selected {dockings}"
            )

        # Get the selector plugin to use
        docking = cast(DockingBase, self.plugins[dockings[0]])

        return docking.run(**kwargs)

    # Finding PDBs for ligands in a folder
    def find_pdb_ligands(self, current_generation_pdb_dir):
        """
        This finds all the pdb files of ligands in a directory

        Inputs:
        :param str current_generation_pdb_dir: the dir path which contains the
            pdb files of ligands to be converted

        Returns:
        :returns: list pdbs_in_folder: a list of all PDB's in the dir
        """

        # TODO: Doesn't seem like the right place for this function.

        # make list of every pdb in the current generations pdb folder
        return list(glob.glob(f"{current_generation_pdb_dir}*.pdb"))

    def run_ligand_handling_for_docking(
        self, pdb_file, file_conversion_class_object: ParentPDBQTConverter
    ) -> Optional[str]:
        """
        this function converts the ligands from PDB to PDBQT format. Returns
        NONE if it worked and the name if it failed to convert.

        Inputs:
        :param str pdb_file: the pdb file of a ligand to format, dock and
            score

        Returns:
        :returns: str smile_name: name of smiles if it failed to dock returns
            None if it docked properly
        """

        # TODO: Doesn't seem like the right place for this.

        assert (
            file_conversion_class_object is not None
        ), "file_conversion_class_object must be passed to VinaDocking"

        # convert ligands to pdbqt format
        # log("\nConverting ligand PDB files to PDBQT format...")
        (
            did_it_convert,
            smile_name,
        ) = file_conversion_class_object.convert_ligand_pdb_file_to_pdbqt(pdb_file)

        if not did_it_convert:
            # conversion failed
            return smile_name

        # Conversion pass. Return None
        # only return failed smile_names which will be handled later
        return None

    # Find ligands which converted to PDBQT
    def find_converted_ligands(self, current_generation_pdb_dir):
        """
        This finds all the pdbqt files of ligands in a directory

        Inputs:
        :param str current_generation_pdb_dir: the dir path which contains the
            pdbqt files of ligands to be docked

        Returns:
        :returns: list pdbqts_in_folder: a list of all PDBqt's in the dir
        """

        # TODO: Not the right place for this, and also specific to vina-like
        # ones.

        # make list of every pdbqt in the current generations pdb folder
        return list(glob.glob(f"{current_generation_pdb_dir}*.pdbqt"))

    def rank_and_save_output_smi(
        self,
        params: dict,
        current_generation_dir: str,
        current_gen_int: int,
        smiles_file: str,
        deleted_smiles_names_list,
    ) -> str:
        """
        Given a folder with PDBQT's, rank all the SMILES based on docking
        score (High to low). Then format it into a .smi file. Then save the
        file.

        Inputs:
        :param dict params: params needs to be threaded here because it has the
            paralizer object which is needed within Scoring.run_scoring_common
        :param str current_generation_dir: path of directory of current
            generation
        :param int current_gen_int: the interger of the current generation
            indexed to zero
        :param str smiles_file:  File path for the file with the ligands for
            the generation which will be a .smi file
        :param list deleted_smiles_names_list: list of SMILES which may have
            failed the conversion process

        Returns:
        :returns: str output_ranked_smile_file: the path of the output ranked
            .smi file
        """

        # TODO: Not the right place for this.

        # Get directory string of PDB files for Ligands
        folder_with_pdbqts = f"{current_generation_dir}PDBs{os.sep}"

        # Run any compatible Scoring Function
        postDockedCompoundInfos = Scoring.run_scoring_common(
            params, smiles_file, folder_with_pdbqts
        )

        # Before ranking these we need to handle Pass-Through ligands from the
        # last generation If it's current_gen_int==1 or if
        # params['redock_elite_from_previous_gen'] is True -Both of these states
        # dock all ligands from the last generation so all of the pass-through
        # lig are already in the PDB's folder thus they should be accounted
        # for in smiles_list If params['redock_elite_from_previous_gen'] is False
        # and current_gen_int != 1 - We need to append the scores form the
        # last gen to smiles_list

        # Only add these when we haven't already redocked the ligand
        if (
            self.params["redock_elite_from_previous_gen"] is False
            and current_gen_int != 0
        ):
            # Go to previous generation folder
            prev_gen_num = str(current_gen_int - 1)
            run_folder = self.params["output_directory"]
            previous_gen_folder = f"{run_folder}generation_{prev_gen_num}{os.sep}"
            ranked_smi_file_prev_gen = (
                f"{previous_gen_folder}generation_{prev_gen_num}_ranked.smi"
            )

            # Also check sometimes Generation 1 won't have a previous
            # generation to do this with and sometimes it will
            if (
                current_gen_int != 1
                or os.path.exists(ranked_smi_file_prev_gen) is not False
            ):
                # Shouldn't happen but to be safe.
                self._process_ligand_scores_from_prev_gen(
                    ranked_smi_file_prev_gen,
                    current_generation_dir,
                    current_gen_int,
                    postDockedCompoundInfos,
                )
        # Output format of the .smi file will be: SMILES    Full_lig_name
        # shorthandname   ...AnyCustominfo... Fitness_metric  diversity
        # Normally the docking score is the fitness metric but if we use a
        # Custom metric than dock score gets moved to index -3 and the new
        # fitness metric gets -2

        # sort list by the affinity of each sublist (which is the last index
        # of sublist)
        postDockedCompoundInfos.sort(key=lambda x: x.score, reverse=False)

        # score the diversity of each ligand compared to the rest of the
        # ligands in the group this adds on a float in the last column for the
        # sum of pairwise comparisons the lower the diversity score the more
        # unique a molecule is from the other mols in the same generation
        postDockedCompoundInfos = Ranking.score_and_calc_diversity_scores(
            postDockedCompoundInfos
        )

        # name for the output file
        output_ranked_smile_file = smiles_file.replace(".smi", "") + "_ranked.smi"

        # save to a new output smiles file. ie. save to ranked_smiles_file

        with open(output_ranked_smile_file, "w") as output:
            for postDockedCompoundInfo in postDockedCompoundInfos:
                output_line = "\t".join(postDockedCompoundInfo.to_list()) + "\n"
                output.write(output_line)

        return output_ranked_smile_file

    def _process_ligand_scores_from_prev_gen(
        self,
        ranked_smi_file_prev_gen: str,
        current_generation_dir: str,
        current_gen_int: int,
        smiles_list: list[PostDockedCompoundInfo],
    ):
        # Note that this modifies the smiles_list in place

        # TODO: Not the right place for this...

        # CHECKED: smiles_list is of type List[PostDockedCompoundInfo] here.

        print("Getting ligand scores from the previous generation")

        # Shouldn't happen but to be safe.
        if os.path.exists(ranked_smi_file_prev_gen) is False:
            raise Exception(
                "Previous generation ranked .smi file does not exist. "
                + "Check if output folder has been moved"
            )

        # Get the data for all ligands from previous generation ranked
        # file
        prev_gen_data_list = Ranking.get_usable_format(ranked_smi_file_prev_gen)
        # CHECKED: prev_gen_data_list of type List[PreDockedCompoundInfo] here.

        # Get the list of pass through ligands
        current_gen_pass_through_smi = (
            current_generation_dir
            + f"SeedFolder{os.sep}Chosen_Elite_To_advance_Gen_{current_gen_int}.smi"
        )
        pass_through_list = Ranking.get_usable_format(current_gen_pass_through_smi)
        # CHECKED: pass_through_list is of type List[PreDockedCompoundInfo] here.

        # Convert lists to searchable Dictionaries.
        prev_gen_data_dict = Ranking.convert_usable_list_to_lig_dict(prev_gen_data_list)
        # CHECKED: prev_gen_data_dict is of type Dict[str, PreDockedCompoundInfo] here.

        assert prev_gen_data_dict is not None, "prev_gen_data_dict is None"

        pass_through_data: List[PostDockedCompoundInfo] = []
        for lig in pass_through_list:
            # CHECKED: lig is of type PreDockedCompoundInfo here.

            lig_data = prev_gen_data_dict[str(lig.smiles + lig.name)]
            # CHECKED: lig_data is of type PreDockedCompoundInfo here.

            # NOTE: Here it must be converted to a PostDockedCompoundInfo
            assert (
                lig_data.previous_docking_score is not None
            ), "lig_data.previous_docking_score is None"

            # TODO: Nervous that additional_info = "". Not sure what to put there.
            lig_info_remove_diversity_info = PostDockedCompoundInfo(
                smiles=lig.smiles,
                id=lig.name,
                short_id=lig.name,
                additional_info="",
                score=lig_data.previous_docking_score,
                diversity_score=None,
            )
            pass_through_data.append(lig_info_remove_diversity_info)

        smiles_list.extend(pass_through_data)
