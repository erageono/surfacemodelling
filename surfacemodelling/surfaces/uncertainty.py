import numpy as np


class UncertaintyVisualization:
    @staticmethod
    def distance(surface, originalpoints, threshold=False, color=False, in_place=True):
        distances = []
        for surfpt in surface.points:
            localdistances = []
            for origpt in originalpoints:
                sfpt_to_orpt = (
                    (origpt[0] - surfpt[0]) ** 2
                    + (origpt[1] - surfpt[1]) ** 2
                    + (origpt[2] - surfpt[2]) ** 2
                ) ** 0.5
                localdistances.append(sfpt_to_orpt)
            distances.append(min(localdistances))

        if in_place:
            surface["scalars"] = np.asarray(distances)
            if threshold is not False:
                surface.clip_scalar(inplace=True, value=threshold)
            return surface
        else:
            copy = surface
            copy["scalars"] = np.asarray(distances)
            if threshold is not False:
                threshed = copy.clip_scalar([0, threshold])
                copy = threshed
            return copy

    @staticmethod
    def uk_stddev(surface,stddev,threshold=False, in_place = True):
        stddev = np.asarray(stddev)
        if in_place:
            surface["scalars"] = stddev
            if threshold is not False:
                surface.clip_scalar(inplace=True, value=threshold)
            return surface
        else:
            copy = surface
            copy["scalars"] = stddev
            if threshold is not False:
                threshed = copy.clip_scalar([0, threshold])
                copy = threshed
            return copy
        

