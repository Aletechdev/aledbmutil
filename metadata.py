import csv


def get_condition_val_str(metadata_row_list):
    l = list(metadata_row_list.copy())
    del l[0]
    condition_str = " ".join(l)
    condition_str = condition_str.replace('"', '')
    condition_str = condition_str.replace(',', ' ')
    return condition_str


def get_condition_val_dict(metadata_file_path):
    condition_val_dict = {}
    with open(metadata_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in reader:
            condition_val_str = get_condition_val_str(row)
#            if '?' not in condition_val_str and 'none' not in condition_val_str:
            condition_val_dict[row[0].lower()] = condition_val_str  # Some labels don't make sense if lowercase, like O2.
                # condition_val_dict[row[0].lower()] = condition_val_str.lower()
    return condition_val_dict


def get_condition_field_val_set(exp_metadata_dict):
    condition_set = set()
    for key_val in exp_metadata_dict.items():
        condition_str = " ".join(key_val)
        condition_str = condition_str.replace('"', '')
        condition_str = condition_str.replace(',', ' ')
        condition_set.add(condition_str)
    return condition_set


import os

def get_all_exp_cond_d(metadata_path):
    all_exp_cond_d = dict()
    
    all_exp_cond_d = {}
    for root, dirs, files in os.walk(metadata_path):
        for name in files:
            relative_file_path_str = str(root)+'/'+name
            d = get_condition_val_dict(relative_file_path_str)
            all_exp_cond_d[d["project"]] = d.copy()
            
    # Removing poor annotation
    no_volume_str = "(0)"
    for m_d in all_exp_cond_d.values():
        for key, value in m_d.items():
            if no_volume_str in value:
                # Works because I'm changing reference.
                m_d[key] = value.replace(no_volume_str, "")
                
    # Adjusting internals of all_exp_cond_d
    # TODO: would probably be nicer code if I instead said which of these to keep.
    metadata_field_remove_l = [
        "Flask-number",
        "ALE-number",
        "Isolate-number",
        "technical-replicate-number",
        "biological-replicates",
        "creator-email",
        "data-type",
        "technical-replicates",
        "serial-number",
        "read-files",
        "link-to-reference-sequence",
        "archive-link",
        # Global removals that could rather be conditional.
        "cultivation-details",
        "pre-culture-details",
        "environment",
        "antibody",
        "antibiotic",
        "isolate-type",
        "experiment-details",
        "sample-time",
        "run-date",
        "project",
        "creator",
        "machine",
        "read-length",
        "read-type",
        "library-prep-kit",
        "library-prep-kit-cycles",
        "library-prep-kit-manufacturer",
        'Link-to-reference-sequence',
        "genome-reference-file",
        "electron-acceptor"
    ]
    for exp_metadata_d in all_exp_cond_d.values():
        for metadata_field_remove_str in metadata_field_remove_l:
            exp_metadata_d.pop(metadata_field_remove_str.lower(), None)
            exp_metadata_d.pop(metadata_field_remove_str, None)

    # improving annotations temperature metadata annotations
    for exp_metadata_d in all_exp_cond_d.values():
        exp_metadata_d["temperature"] += " celsius"
        #     exp_metadata_d["strain"] = str(exp_metadata_d["taxonomy-id"]+' '+exp_metadata_d["strain-description"])
        #     del exp_metadata_d["taxonomy-id"]
        #     del exp_metadata_d["strain-description"]

    return all_exp_cond_d
