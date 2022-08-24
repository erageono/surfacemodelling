import csv
import glob
import os
import geojson

class Files:
    @staticmethod
    def kof_points(path):
        """Les *.kof fil med punkter

        Les alle *.kof filer fra en spesifert mappe
        Read all *.kof files from a specified folder

        Args:
            path (str): Linke til data mappe

        Returns:
            pointslist (list): Liste med format [Navn (str), x (float), y (float), z (float), ]
        """
        pointslist = []
        for filename in glob.glob(os.path.join(path, "*.kof")):
            with open(filename, newline="") as infile:
                inreader = csv.reader(infile, delimiter=" ",)
                for row in inreader:
                    cleanrow = list(filter(None, row))
                    if cleanrow[0] == "05":
                        pointslist.append(
                            [
                                float(cleanrow[4]),
                                float(cleanrow[3]),
                                float(cleanrow[5]),
                            ]
                        )
        return pointslist

    @staticmethod
    def kof_file(filename):
        """Les *.kof fil med punkter

        Args:
            path (str): Linke til .kof fil

        Returns:
            pointslist (list): Liste med format [x (float), y (float), z (float), ]
        """
        pointslist = []
        with open(filename, newline="") as infile:
            inreader = csv.reader(infile, delimiter=" ",)
            for row in inreader:
                cleanrow = list(filter(None, row))
                if cleanrow[0] == "05":
                    pointslist.append(
                        [
                            float(cleanrow[4]),
                            float(cleanrow[3]),
                            float(cleanrow[5]),
                        ]
                    )
        return pointslist

    @staticmethod
    def geojson(path):
        """Les *.geojson fil med punkter og instruksjoner/andre egenskaper om overflatter

        Args:
            path (str): Linke til fil

        Returns:
            instructions (list): Dictionary med format f.eks:
                                 {{'Berg': {'color': '#000000', 'type': 'Layer','points': [...List of source points...] 'filepath': False, 'interpretation': {'outputCRS': 'UTM33_NN2000', 'method': 'gempy.decimated_smoothed', 'priority': 1, 'buffer': [100, 100, 100], 'uncertainty': {'method': 'distance', 'threshold': 1000}}},...
        """  # noqa E501
        instructions = {}
        with open(path) as f:
            gj = geojson.load(f)
        features = gj["features"]
        for feature in features:
            featurename = feature["properties"]["identity"]["layerName"]
            instructions[featurename] = {
                'color': feature["properties"]["identity"]["layerColor"],
                'type': 'Layer',
                'points': feature["geometry"]["coordinates"],
            }

            if (
                len(feature["properties"]["meshfilepath"]) > 3
                and type(feature["properties"]["meshfilepath"]) == str
            ):
                instructions[featurename]["filepath"] = feature["properties"]["meshfilepath"]
            else:
                instructions[featurename]["filepath"] = False
            instructions[featurename]["interpretation"] = {}
            instructions[featurename]["interpretation"]["outputCRS"] = feature["properties"][
                "outputCRS"
            ]
            instructions[featurename]["interpretation"]["method"] = feature["properties"][
                "interpolation"
            ]["interpolationMethod"]
            instructions[featurename]["interpretation"]["priority"] = feature["properties"][
                "interpolation"
            ]["layerPriority"]
            instructions[featurename]["interpretation"]["buffer"] = [
                feature["properties"]["interpolation"]["interpolationBuffer"]["x"],
                feature["properties"]["interpolation"]["interpolationBuffer"]["y"],
                feature["properties"]["interpolation"]["interpolationBuffer"]["z"],
            ]
            if feature["properties"]["uncertainty"] != {}:
                instructions[featurename]["interpretation"]["uncertainty"] = {}
                instructions[featurename]["interpretation"]["uncertainty"]["method"] = feature[
                    "properties"
                ]["uncertainty"]["uncertaintyMethod"]
                instructions[featurename]["interpretation"]["uncertainty"]["threshold"] = feature[
                    "properties"
                ]["uncertainty"]["uncertaintyDisplayThreshold"]
            else:
                instructions[featurename]["interpretation"]["uncertainty"] = False
        return instructions

    @staticmethod
    def default_cosite_kof(path):
        """Les flere *.kof fil med punkter hvor en fil har punkter fra en lag

        Les alle *.kof filer fra en spesifert mappe
        Read all *.kof files from a specified folder

        Args:
            path (str): Lenke til data mappen

        Returns:
            pointslist (list): Liste med format [Navn (str), x (float), y (float), z (float), ]
        """
        
        layers = {}
        for filename in glob.glob(os.path.join(path, "*.kof")):
            layername = filename.strip("tolkning.kof")
            layername = layername.split("-")[-2]
            print(layername)
            layers[layername] = []
            with open(filename, newline="") as infile:
                inreader = csv.reader(
                    infile,
                    delimiter=" ",
                )
                for row in inreader:
                    cleanrow = list(filter(None, row))
                    if cleanrow[0] == "05":
                        layers[layername].append(
                            [
                                float(cleanrow[4]),
                                float(cleanrow[3]),
                                float(cleanrow[5]),
                            ]
                        )
        return layers