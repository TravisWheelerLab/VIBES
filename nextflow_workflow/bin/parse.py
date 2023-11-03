#!/usr/bin/env python3
import argparse
import glob
import json
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Produce SODA-based visualizations from the output of the VIBES pipeline"
    )

    parser.add_argument(
        type=str,
        help="The top-level directory of VIBES output files",
        dest="vibes_output_dir",
        metavar="<dir>",
    )

    parser.add_argument(
        "-b",
        type=str,
        help="Path to the VIBES-SODA JavaScript bundle",
        metavar="<bundle.js>",
        dest="bundle",
        default="./vibes-soda.js",
    )

    parser.add_argument(
        "-t",
        type=str,
        help="Path to the template HTML file",
        metavar="<template.html>",
        dest="template",
        default="./template.html",
    )

    parser.add_argument(
        "-o",
        type=str,
        help="Path to the output directory",
        metavar="<dir>",
        dest="outdir",
        default="./viz",
    )

    args = parser.parse_args()

    return args


def quote_str(string: str) -> str:
    return f'"{string}"'


def name_from_path(path: str) -> str:
    path_tokens = path.split("/")
    name_tokens = path_tokens[-1].split(".")
    # just in case the file name has any "." in it
    name = ".".join(name_tokens[:-1])

    return name


class Integration:
    def __init__(self, record: str):
        tokens = record.split("\t")

        start = int(tokens[10])
        end = int(tokens[11])

        if start > end:
            tmp = start
            start = end
            end = tmp

        self.query_name = tokens[0]
        self.accession = tokens[2]
        self.evalue = float(tokens[3])
        self.full_length = bool(tokens[4])
        self.query_start = int(tokens[5])
        self.query_end = int(tokens[6])
        self.query_length = int(tokens[7])
        self.target_name = tokens[9]
        self.target_start = start
        self.target_end = end
        self.target_length = int(tokens[12])
        self.strand = tokens[13]

    def to_string(self):
        return "{},{},{},{},{},{},{},{}".format(
            self.target_start,
            self.target_end,
            self.query_start,
            self.query_end,
            self.strand,
            self.evalue,
            self.query_name,
            self.accession,
        )


class BacterialGene:
    def __init__(self, record: str):
        # gnl|Prokka|NDHMLNHN_1
        # Prodigal:002006
        # CDS
        # 1
        # 1371
        # .
        # +
        # 0
        #   ID=NDHMLNHN_00001;
        #   Parent=NDHMLNHN_00001_gene;
        #   eC_number=3.6.-.-;
        #   Name=mnmE;
        #   db_xref=COG:COG0486;
        #   gene=mnmE;
        #   inference=ab initio prediction:Prodigal:002006,similar to AA sequence:UniProtKB:P25522;
        #   locus_tag=NDHMLNHN_00001;
        #   product=tRNA modification GTPase MnmE;
        #   protein_id=gnl|Prokka|NDHMLNHN_00001

        tokens = record.split("\t")

        meta_dict = {}
        meta_tokens = tokens[8].split(";")
        for meta_token in meta_tokens:
            [key, val] = meta_token.split("=")
            meta_dict[key] = val

        start = int(tokens[3])
        end = int(tokens[4])

        if start > end:
            tmp = start
            start = end
            end = tmp

        self.target_name = tokens[0]
        self.target_start = start
        self.target_end = end
        self.strand = tokens[6]

        if "locus_tag" in meta_dict:
            self.id = meta_dict["locus_tag"]
        else:
            self.id = ""

        try:
            self.score = float(tokens[7])
        except ValueError:
            self.score = 0.0

        if "product" in meta_dict:
            self.product = meta_dict["product"]
        else:
            self.product = ""

        if "gene" in meta_dict:
            self.query_name = meta_dict["gene"]
        elif "Name" in meta_dict:
            self.query_name = meta_dict["Name"]
        else:
            self.query_name = self.id

    def to_string(self):
        return "{},{},{},{},{},{},{}".format(
            self.target_start,
            self.target_end,
            self.strand,
            self.score,
            self.query_name,
            self.id,
            self.product,
        )


class ViralGene:
    def __init__(self, record: str):
        tokens = record.split("\t")

        start = int(tokens[10])
        end = int(tokens[11])

        if start > end:
            tmp = start
            start = end
            end = tmp

        self.query_name = tokens[0]
        self.description = tokens[1]
        self.accession = tokens[2]
        self.evalue = float(tokens[3])
        self.full_length = bool(tokens[4])
        self.query_start = int(tokens[5])
        self.query_end = int(tokens[6])
        self.query_length = int(tokens[7])
        self.target_name = tokens[9]
        self.target_start = start
        self.target_end = end
        self.target_length = int(tokens[12])
        self.strand = tokens[13]

    def to_string(self):
        return "{},{},{},{},{},{},{},{},{},{}".format(
            self.target_start,
            self.target_end,
            self.query_start,
            self.query_end,
            self.query_length,
            self.strand,
            self.evalue,
            self.query_name,
            self.accession,
            self.description,
        )


class Sequence:
    def __init__(self, name: str, length: int):
        self.name = name
        self.length = length


class Occurrence:
    def __init__(self, integrations: [Integration]):
        self.name = integrations[0].query_name
        self.counts = [0 for _ in range(integrations[0].query_length)]
        for integration in integrations:
            for i in range(integration.query_start, integration.query_end):
                self.counts[i] += 1


def parse_integration_tsv(path: str) -> [Integration]:
    integrations = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#"):
                continue

            integrations.append(Integration(line))

    return integrations


def parse_viral_gene_tsv(path: str) -> [ViralGene]:
    genes = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#"):
                continue

            genes.append(ViralGene(line))

    return genes


def parse_bacterial_gene_gff3(path: str) -> ([ViralGene], [Sequence]):
    genes = []
    sequences = []
    with open(path, "r") as f:
        lines = f.readlines()

        header_lines = [l for l in lines if l.startswith("##sequence-region")]

        if len(header_lines) == 0:
            print(f"no sequence-region header lines found in gff file: {path}")
            exit()

        for line in header_lines:
            tokens = line.split(" ")
            sequences.append(Sequence(tokens[1], int(tokens[3])))

        record_lines = [l for l in lines if not l.startswith("#")]

        for line in record_lines:
            if line.startswith(">"):
                # we've reached the end of
                # the records at this point
                break

            genes.append(BacterialGene(line))

    return (genes, sequences)


def parse_records(args):
    integrations_paths = glob.glob(
        f"./{args.vibes_output_dir}/tsv/bacterial_integrations/*.tsv"
    )

    viral_gene_paths = glob.glob(
        f"./{args.vibes_output_dir}/tsv/viral_gene_annotations/*.tsv"
    )

    bacterial_gene_paths = glob.glob(f"./{args.vibes_output_dir}/gff/*.gff")

    integrations_paths.sort()
    viral_gene_paths.sort()
    bacterial_gene_paths.sort()

    # we pull the bacteria names from the file names,
    # since all of the records refer to the individual
    # sequences in the bacterial genomes
    bacteria_names: [str] = []

    ##
    all_integrations: [Integration] = []
    for path in integrations_paths:
        bacteria_names.append(name_from_path(path))
        all_integrations += parse_integration_tsv(path)

    ##
    all_viral_genes: [ViralGene] = []
    for path in viral_gene_paths:
        all_viral_genes += parse_viral_gene_tsv(path)

    ##
    all_bacterial_genes: [BacterialGene] = []
    sequence_map: {str: [Sequence]} = {}
    for path in bacterial_gene_paths:
        bacteria_name = name_from_path(path)
        bacteria_names.append(bacteria_name)

        (genes, seqs) = parse_bacterial_gene_gff3(path)
        all_bacterial_genes += genes
        sequence_map[bacteria_name] = seqs

    # we grabbed bacteria names from both integration
    # and gene files, so we filter for unique names
    bacteria_names = list(set(bacteria_names))

    # grab all virus names by taking the unique
    # query names amongst all integrations
    virus_names = list(set([i.query_name for i in all_integrations]))

    # once we have all the data parsed,
    # we can build the occurrences
    all_occurrences: [Occurrence] = []
    for name in virus_names:
        filtered = [
            integration
            for integration in all_integrations
            if integration.query_name == name
        ]

        all_occurrences.append(Occurrence(filtered))

    data_list = []

    for bacteria_name in bacteria_names:
        sequences = []
        sequence_names = []
        for s in [s for s in sequence_map[bacteria_name]]:
            sequence_names.append(s.name)
            sequences.append(
                {
                    "sequenceName": s.name,
                    "sequenceLength": s.length,
                    "integrations": [
                        i.to_string()
                        for i in all_integrations
                        if i.target_name == s.name
                    ],
                    "genes": [
                        g.to_string()
                        for g in all_bacterial_genes
                        if g.target_name == s.name
                    ],
                }
            )

        virus_names = list(
            set(
                [
                    i.query_name
                    for i in all_integrations
                    if i.target_name in sequence_names
                ]
            )
        )

        occurrences = []
        for virus in virus_names:
            occ = [o for o in all_occurrences if o.name == virus]
            assert len(occ) == 1

            occurrences.append(
                {
                    "virusName": virus,
                    "counts": occ[0].counts,
                    "genes": [
                        g.to_string() for g in all_viral_genes if g.target_name == virus
                    ],
                }
            )

        data = {
            "bacteriaName": bacteria_name,
            "virusData": occurrences,
            "bacteriaData": sequences,
        }

        data_list.append(data)

    return data_list


def write_data(args, data_list):
    # html stuff
    vibes_soda_bundle = open(f"{args.bundle}").read()
    template_html = open(f"{args.template}").read()
    viz_html = template_html.replace("VIBES_SODA_TARGET", vibes_soda_bundle)

    # root output directory
    os.makedirs(f"{args.outdir}", exist_ok=True)

    # these are the names of the input bacterial genome fasta files,
    # i.e. this is usually something like "Pseudomonas_blah_blah_blah_number"
    bacteria_names = [l["bacteriaName"] for l in data_list]

    # sentinel file
    with open(f"{args.outdir}/bacteria.js", "w") as out:
        out.write(
            "bacteriaNames = [{}];".format(
                ",".join([quote_str(n) for n in bacteria_names])
            )
        ),

    for data in data_list:
        bacteria_name = data["bacteriaName"]
        json_str = json.dumps(data)
        data_str = f"let data = {json_str};"
        with open(f"{args.outdir}/{bacteria_name}.html", "w") as out:
            out.write(viz_html.replace("VIBES_DATA_TARGET", data_str))

    pass


def main():
    args = parse_args()
    data_list = parse_records(args)
    write_data(args, data_list)


main()
