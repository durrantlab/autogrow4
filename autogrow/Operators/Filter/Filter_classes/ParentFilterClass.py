class ParentFilter(object):
    """
    This is a scipt containing all of the filters for drug likeliness

    Filters for orally bio-available drugs:
        1) Lipinski

    Filters for for lead-likeness:
        1) Ghose
        2) Mozziconacci

    Filters for CNS/Blood Brain Barrier Permeable:
        1) VandeWaterbeemd

    False-Positive/Metabolite substructure searches:
        1) PAINS_Filter
        2) NIH_Filter
        3) BRENK_Filter
    """
    def get_name(self):
        """
        Returns the current class name.    
        Returns:
        :returns: str self.__class__.__name__: the current class name.

        """
        return self.__class__.__name__
    #
    def run_filter(self, input_string):
        """
        run_filter is needs to be implimented in each class.
        Inputs:
        :param str input_string:  A string to raise an exception
        """
        raise NotImplementedError("run_filter() not implemented")
    #
    