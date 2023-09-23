#!/bin/python3
def sample_grabber(input):
    ''' This function takes as input the project.list and returns a list of samples within the given project. '''
    prefix_list = []
    with open(input, 'r') as prefix_file:
        header = prefix_file.readline()
        for line in prefix_file:
            columns = line.strip().split('\t')
            prefix_list.append(columns[0])
    return(prefix_list)

def compare_grabber(DE_dir, project):
    ''' This function finds and retuns all the comparisons that the user has outlined they want to run under project_dir/DE_analysis '''
    compare = []
    suffix = "_gtf.list"
    prefix = project + "_"
    files = os.listdir(DE_dir)
    for file in files:
        if file.endswith(suffix):
            new_name = file.replace(suffix, "")
            new_name = new_name.replace(prefix, "")
            compare.append(new_name)
    for cond in compare:
        compare_dir = os.path.join(DE_dir, cond)
        if not os.path.exists(compare_dir):
            os.makedirs(compare_dir)
    return(compare)

def make_executable(path):
    ''' This function takes an input a script and makes it executable '''
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2    # copy R bits to X
    os.chmod(path, mode)

def DE_comparisons(sample_list, project_dir, project, DE_dir):
    '''  '''
    lists_dict = {}
    #Collect all unique identifiers and create key:value pairs for them
    with open(sample_list, 'r') as samples:
        #Skip header
        next(samples)
        for line in samples:
            col = line.strip().split('\t')
            #Get column values
            sample_id = col[0]
            identifier = col[1].split(',')
            for ID in identifier:
                #Create a unique list name based off identified and type_value
                list_name = ID
                #Check if the list is new
                if list_name not in lists_dict:
                    lists_dict[list_name] = []
                #Add matched sample_id to correct list
                lists_dict[list_name].append(sample_id)

    #Based off all unique key pairs create project_comparisons.list for downstream DE
    comparison_file = os.path.join(project_dir, project + "_comparisons.list")
    with open(comparison_file, 'w') as comparisons:
        comparisons.write("Subproject" + '\t')
        comparisons.write("Treated" + '\t')
        comparisons.write("Control" + '\n')
        keys = list(lists_dict.keys())
        for i in range(len(keys)):
            key1 = keys[i]
            for j in range(i+1, len(keys)):
                key2 = keys[j]
                comparisons.write(f"{key1}_v_{key2}\t{key1}\t{key2}\n")


    #Read in newly created project_comparisons.list and based off the key:values create the gtf list files to be passsed to downstream steps
    with open(comparison_file, 'r') as comparisons:
        #Skip header
        next(comparisons)
        for line in comparisons:
            treatment_ID = []
            control_ID = []
            col = line.strip().split('\t')
            compare = col[0]
            treated = col[1]
            control = col[2]
            compare_path = os.path.join(DE_dir, compare)
            if not os.path.exists(compare_path):
                print(f"Creating subfolder: {compare}")
                os.makedirs(compare_path)
            for key,value in lists_dict.items():
                if key == treated:
                    treatment_ID.extend(value)
                elif key == control:
                    control_ID.extend(value)
            compare_file_gtf = compare + "_gtf.list"
            compare_file_geneout = compare + "_geneout.list"
            compare_file_gtf_path = os.path.join(DE_dir, compare, compare_file_gtf)
            compare_file_geneout_path = os.path.join(DE_dir, compare, compare_file_geneout)
            with open(compare_file_geneout_path, 'w') as geneout_out, open(compare_file_gtf_path, 'w') as gtf_out:
                for sample in treatment_ID:
                    sample = str(sample)
                    gtf_path = project_dir + "/" + sample + "/circular/output/" + sample + ".gtf"
                    geneout_path = project_dir + "/" + sample + "/circular/output/gene/" + sample + "_out.gtf"
                    gtf_out.write(sample + '\t')
                    gtf_out.write(gtf_path + '\t')
                    gtf_out.write("T" + '\n')
                    geneout_out.write(sample + '\t')
                    geneout_out.write(geneout_path + '\n')
                for sample in control_ID:
                    sample = str(sample)
                    gtf_path = project_dir + "/" + sample + "/circular/output/" + sample + ".gtf"
                    geneout_path = project_dir + "/" + sample + "/circular/output/gene/" + sample + "_out.gtf"
                    gtf_out.write(sample + '\t')
                    gtf_out.write(gtf_path + '\t')
                    gtf_out.write("C" + '\n')
                    geneout_out.write(sample + '\t')
                    geneout_out.write(geneout_path + '\n')

def compare_grabber(DE_dir, project):
    ''' This function finds and retuns all the comparisons that the user has outlined they want to run under project_dir/DE_analysis '''
    compare = []
    for folder_name in os.listdir(DE_dir):
        folder_path = os.path.join(DE_dir, folder_name)
        if os.path.isdir(folder_path):
            compare.append(folder_name)
    return(compare)

def initilize_DE(pipeline_dir, DE_dir, project, compare_list):
    ''' This function when called runs all 3 scripts that make up the DE analysis '''
    pipeline = pipeline_dir.rstrip("/")
    script1 = os.path.join(pipeline, "DE_analysis_step1_rep.sh")
    script2 = os.path.join(pipeline, "DE_analysis_step1_norep.sh")
    script3 = os.path.join(pipeline, "DE_analysis_step2.py")
    script4 = os.path.join(pipeline, "DE_analysis_step3.sh")
    for compare in compare_list:
        T_counter = 0
        C_counter = 0
        check_file = compare + "_gtf.list"
        check_path = os.path.join(DE_dir, compare, check_file)
        compare_dir = os.path.join(DE_dir, compare)
        with open(check_path, 'r') as mode_check:
            for line in mode_check:
                items = line.strip().split()
                if items[2] == "T":
                    T_counter += 1
                elif items[2] == "C":
                    C_counter += 1
            if (T_counter < 1) | (C_counter < 1):
                print(f"Invalid configuation for analysis DE, please verify {compare}")
                run_mode = "invalid"
            elif (T_counter == 1) & (C_counter == 1):
                print(f"Single sample comparison DE analysis will be ran on {compare}")
                run_mode = "single_DE"
            elif (T_counter > 1) | (C_counter > 1):
                print(f"Multiple replicate comparison DE analysis will be ran on {compare}")
                run_mode = "replicate_DE"
        if run_mode == "single_DE":
            #Run CIRI_DE
            with open(check_path, 'r') as mode_check:
                for line in mode_check:
                    items = line.strip().split()
                    if items[2] == "T":
                        T_sample = items[1]
                    elif items[2] == "C":
                        C_sample = items[1]
            CIRI_DE_output = os.path.join(compare_dir, compare + "_DE.csv")
            command = "bash " + script2 + " " + C_sample + " " + T_sample + " " + CIRI_DE_output
            process = subprocess.run(command, shell=True)

        elif run_mode == "replicate_DE":
            #Run prep_CIRIquant
            prefix = project + "_" + compare
            output = compare_dir + "/"
            command = "bash " +  script1 + " " + check_path + " " + prefix + " " + output
            process = subprocess.run(command, shell=True)
            #Run prepDE.py
            list = DE_dir + "/" + project + "_" + compare + "_geneout.list"
            gene = compare_dir + "/" + project + "_" + compare + "_gene_count_matrix.csv"
            transcript = compare_dir + "/" + project + "_" + compare + "_transcript_count_matrix.csv"
            command = "python2 " + script3 + " -i " + list + " -g " + gene + " -t " + transcript
            process = subprocess.run(command, shell=True)
            #Run CIRI_DE_replicate
            pre = compare_dir + "/" + project + "_" + compare
            command = "bash " + script4 + " " + pre
            process = subprocess.run(command, shell=True)



