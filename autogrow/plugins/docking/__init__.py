from abc import abstractmethod
import glob
import os
import random
from typing import Any, Dict, List, Optional, Tuple, Union, cast
from autogrow.config.argparser import ArgumentVars
from autogrow.plugins.plugin_base import PluginBase
from autogrow.plugins.plugin_manager_base import PluginManagerBase
from autogrow.types import PostDockedCompound, PreDockedCompound, ScoreType
import autogrow.docking.scoring.execute_scoring_mol as Scoring
import autogrow.docking.ranking.ranking_mol as Ranking
from autogrow.utils.logging import log_debug


class DockingBase(PluginBase):
    def run(self, **kwargs) -> List[PostDockedCompound]:
        """Run the plugin with provided arguments."""
        return self.run_docking(predocked_cmpds=kwargs["predocked_cmpds"])

    @abstractmethod
    def run_docking(
        self, predocked_cmpds: List[PreDockedCompound]
    ) -> List[PostDockedCompound]:
        """
        run_docking is needs to be implemented in each class.

        Inputs:
        :param PreDockedCompound predocked_cmpds: A List of PreDockedCompound
            objects.

        Returns:
        :returns: List[PostDockedCompound]: A list of PostDockedCompound
            objects, each containing the score and a docked (posed) SDF file.
        """

        # raise NotImplementedError("run_dock() not implemented")
        pass


class DockingPluginManager(PluginManagerBase):
    def run(self, **kwargs) -> List[PostDockedCompound]:
        """
        Run the plugin with provided arguments.

        Inputs:
        :param dict kwargs: a dictionary of arguments to pass to the plugin

        Returns:
        :returns: List[PostDockedCompound]: A list of PostDockedCompound
            objects, each containing the score and a docked (posed) SDF file.
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

        post_docked_cmpds = docking.run(**kwargs)

        for post_docked_cmpd in post_docked_cmpds:
            log_debug(
                f"Docked molecule {post_docked_cmpd.smiles}. Score: {post_docked_cmpd.docking_score:.2f}"
            )

        return post_docked_cmpds

    def rank_and_save_output_smi(
        self,
        current_generation_dir: str,
        current_gen_int: int,
        smiles_file: str,
        postDockedCompoundInfos: List[PostDockedCompound],
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

        Returns:
        :returns: str output_ranked_smile_file: the path of the output ranked
            .smi file
        """

        # TODO: Not the right place for this.

        # Get directory string of PDB files for Ligands
        # folder_with_pdbqts = f"{current_generation_dir}PDBs{os.sep}"

        # Run any compatible Scoring Function
        # postDockedCompoundInfos = Scoring.run_scoring_common(
        #     params, smiles_file, folder_with_pdbqts
        # )

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
        postDockedCompoundInfos.sort(key=lambda x: x.docking_score, reverse=False)

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
                lst = postDockedCompoundInfo.to_list()
                lst[-1] = os.path.basename(lst[-1])
                output_line = "\t".join(lst) + "\n"
                output.write(output_line)

        return output_ranked_smile_file

    def _process_ligand_scores_from_prev_gen(
        self,
        ranked_smi_file_prev_gen: str,
        current_generation_dir: str,
        current_gen_int: int,
        smiles_list: list[PostDockedCompound],
    ):
        # Note that this modifies the smiles_list in place

        # TODO: Not the right place for this...

        # CHECKED: smiles_list is of type List[PostDockedCompound] here.

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
        # CHECKED: prev_gen_data_list of type List[PreDockedCompound] here.

        # Get the list of pass through ligands
        current_gen_pass_through_smi = (
            current_generation_dir
            + f"SeedFolder{os.sep}Chosen_Elite_To_advance_Gen_{current_gen_int}.smi"
        )
        pass_through_list = Ranking.get_usable_format(current_gen_pass_through_smi)
        # CHECKED: pass_through_list is of type List[PreDockedCompound] here.

        # Convert lists to searchable Dictionaries.
        prev_gen_data_dict = Ranking.convert_usable_list_to_lig_dict(prev_gen_data_list)
        # CHECKED: prev_gen_data_dict is of type Dict[str, PreDockedCompound] here.

        assert prev_gen_data_dict is not None, "prev_gen_data_dict is None"

        pass_through_data: List[PostDockedCompound] = []
        for lig in pass_through_list:
            # CHECKED: lig is of type PreDockedCompound here.

            lig_data = prev_gen_data_dict[str(lig.smiles + lig.name)]
            # CHECKED: lig_data is of type PreDockedCompound here.

            # NOTE: Here it must be converted to a PostDockedCompound
            assert (
                lig_data.previous_docking_score is not None
            ), "lig_data.previous_docking_score is None"

            # TODO: Nervous that additional_info = "". Not sure what to put there.
            lig_info_remove_diversity_info = PostDockedCompound(
                smiles=lig.smiles,
                id=lig.name,
                short_id=lig.name,
                additional_info="",
                docking_score=lig_data.previous_docking_score,
                diversity_score=None,
            )
            pass_through_data.append(lig_info_remove_diversity_info)

        smiles_list.extend(pass_through_data)
