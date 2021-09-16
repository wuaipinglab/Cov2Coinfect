import os
import glob
from collections import defaultdict
from copy import deepcopy

import pandas as pd
import numpy as np
import vcf


def resolve_MNV(temp_records: dict) -> list:
    
    if len(temp_records) == 1:
        (records, ) = temp_records.values()
        return records
    
    temp_records = list(temp_records.values())
    sort_i = list(range(len(temp_records)))
    for records in temp_records:
        for _ in records:
            swap_flag0 = False
            i = 0
            while i < len(records) - 1:
                if records[i]["Frequency"] < records[i + 1]["Frequency"]:
                    records[i], records[i + 1] = records[i + 1], records[i]
                    swap_flag0 = True
                i += 1
            if not swap_flag0:
                break
        j = 0
        while j < len(temp_records) - 1:
            if len(temp_records[j]) < len(temp_records[j + 1]):
                temp_records[j], temp_records[j + 1] = temp_records[j + 1], temp_records[j]
                sort_i[j], sort_i[j + 1] = sort_i[j + 1], sort_i[j]
            j += 1
    
    for _ in temp_records:
        swap_flag1 = False
        k = 0
        while k < len(temp_records) - 1:
            if (temp_records[k][0]["Frequency"] < temp_records[k + 1][0]["Frequency"] and
                len(temp_records[k]) == len(temp_records[k + 1])):
                temp_records[k], temp_records[k + 1] = temp_records[k + 1], temp_records[k]
                sort_i[k], sort_i[k + 1] = sort_i[k + 1], sort_i[k]
                swap_flag1 = True
            k += 1
        if not swap_flag1:
            break
            
    ref_records = temp_records[0]
    res = deepcopy(temp_records[0])
    for records in temp_records[1:]:
        init_compare = 0
        for rec in records:
            min_diff = float('inf')
            compare_start, compare_end = init_compare, init_compare + 1
            ref_freq = 0
            for n, ref_rec in enumerate(ref_records, start=init_compare):
                ref_freq = ref_rec["Frequency"]
                diff = abs(ref_freq - rec["Frequency"])
                if diff < min_diff:
                    min_diff = diff
                    compare_start, compare_end = n, n + 1
            ref_acc_freq = 0
            for n, ref_rec in enumerate(ref_records, start=init_compare):
                ref_acc_freq += ref_rec["Frequency"]
                diff = abs(ref_acc_freq - rec["Frequency"])
                if diff < min_diff:
                    min_diff = diff
                    compare_start, compare_end = init_compare, n + 1
                
            init_compare = compare_end
                    
            for ref_rec_index in range(compare_start, compare_end):
                if type(res[ref_rec_index]["Reference Position"]) is int:
                    res[ref_rec_index]["Reference Position"] = [
                        res[ref_rec_index]["Reference Position"],
                        rec["Reference Position"]
                    ]
                else:
                    res[ref_rec_index]["Reference Position"].append(rec["Reference Position"])

                res[ref_rec_index]["Type"] = "MNV"
                res[ref_rec_index]["Length"] += 1
                res[ref_rec_index]["Reference"] += rec["Reference"]
                res[ref_rec_index]["Allele"] += rec["Allele"]

#                 res[ref_rec_index]["Coverage"] += rec["Coverage"]
#                 res[ref_rec_index]["Frequency"] += rec["Frequency"]
                res[ref_rec_index]["Average quality"] += rec["Average quality"]
    
    for ref_rec in res:
        if (type(ref_rec["Reference Position"]) is int):
            continue
        else:
            sort_i2 = np.argsort(ref_rec["Reference Position"])
            all(i == j for i, j in zip(sort_i, sort_i2))
            ref_rec["Reference Position"] = min(ref_rec["Reference Position"])
            ref_rec["Reference"] = "".join(ref_rec["Reference"][i] for i in sort_i)
            ref_rec["Allele"] = "".join(ref_rec["Allele"][i] for i in sort_i)
    #         ref_rec["Coverage"] /= len(sort_i)
    #         ref_rec["Coverage"] = int(ref_rec["Coverage"])
            ref_rec["Average quality"] /= len(sort_i)
    return res

def vcf2pandas(working_dir: str, vcf_dir: str = None) -> str:
    csv_dir = os.path.join(working_dir, "snp_csv/")
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

    if vcf_dir is None:
        vcf_dir = os.path.join(working_dir, "vcf")

    for vcf_file in glob.glob(os.path.join(vcf_dir, "*.vcf")):
        res = []
        with open(vcf_file) as f:
            prev_pos = -1
            curr_pos = 0
            temp_records = defaultdict(list)
            for record in vcf.Reader(f):
                curr_pos = record.POS
                if "DP4" in record.INFO:
                    total_count = sum(record.INFO["DP4"])
                else:
                    total_count = sum(record.INFO["I16"][:4])
                    
                for allele in record.ALT:
                    reference = record.REF
                    allele_seq = allele.sequence
                    if len(reference) > len(allele_seq):
                        variant_type = "Deletion"
                        reference = reference[len(allele_seq):]
                        variant_len = len(reference)
                        allele_seq = "-"
                        record.POS += 1
                    elif len(record.REF) < len(allele.sequence):
                        variant_type = "Insertion"
                        allele_seq = allele_seq[len(reference):]
                        variant_len = len(allele_seq)
                        reference = "-"
                        record.POS += 1
                    else:
                        variant_len = len(allele_seq)
                        if variant_len == 1:
                            variant_type = "SNV"
                        else:
                            variant_type = "MNV"
                    allele_frequency = record.INFO["AF"]
                    curr_record = {
                        "Reference Position": record.POS,
                        "Type": variant_type,
                        "Length": variant_len,
                        "Reference": reference,
                        "Allele": allele_seq,
                        "Linkage": "",
                        "Count": int(allele_frequency * total_count),
                        "Coverage": int(total_count),
                        "Frequency": allele_frequency * 100,
                        "Forward/reverse balance": record.INFO["SB"],
                        "Average quality": record.QUAL,
                        "Overlapping annotations": "N/A",
                        "Coding region change": "N/A",
                        "Amino acid change": "N/A"
                    }
                    if (curr_pos - prev_pos <= 1 and
                        curr_record["Type"] == "SNV" and
                        prev_record["Type"] == "SNV"):
                        temp_records[curr_record["Reference Position"]].append(curr_record)
                    else:
                        if prev_pos > 0:
                            res.extend(resolve_MNV(temp_records))
                        temp_records = defaultdict(list)
                        temp_records[curr_record["Reference Position"]].append(curr_record)
                prev_record = curr_record
                prev_pos = curr_pos
            if len(temp_records):
                res.extend(resolve_MNV(temp_records))
            res = pd.DataFrame.from_records(res)
            csv_name = os.path.basename(vcf_file)
            csv_name = os.path.splitext(csv_name)[0] + ".csv"
            res.to_csv(os.path.join(csv_dir, csv_name), index=False)
    return csv_dir
    