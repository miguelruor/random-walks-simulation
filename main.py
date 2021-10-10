import classic_simulation
import json
import os
import networkx as nx
import ibm_boto3
from ibm_botocore.client import Config
from ibm_botocore.exceptions import ClientError
from dotenv import load_dotenv

load_dotenv()
# Constants for IBM COS values
COS_ENDPOINT = "https://s3.us-south.cloud-object-storage.appdomain.cloud" 
COS_API_KEY_ID = os.environ['COS_APIKEY'] 
COS_AUTH_ENDPOINT = "https://iam.cloud.ibm.com/identity/token"
COS_SERVICE_CRN = os.environ['COS_SERVICE_CRN']

# Create client connection
cos_cli = ibm_boto3.client("s3",
    ibm_api_key_id=COS_API_KEY_ID,
    ibm_service_instance_id=COS_SERVICE_CRN,
    ibm_auth_endpoint=COS_AUTH_ENDPOINT,
    config=Config(signature_version="oauth"),
    endpoint_url=COS_ENDPOINT
)

def log_done():
    print("DONE!\n")

def log_client_error(e):
    print("CLIENT ERROR: {0}\n".format(e))

def log_error(msg):
    print("UNKNOWN ERROR: {0}\n".format(msg))

# Retrieve a particular item from the bucket
def get_item(bucket_name, item_name):
    print("Retrieving item from bucket: {0}, key: {1}".format(bucket_name, item_name))

    try:
        file = cos_cli.get_object(Bucket=bucket_name, Key=item_name)
        #print("File Contents: {0}".format(file["Body"].read()))
        log_done()

        data = file['Body'].read()

        return data

    except ClientError as be:
        log_client_error(be)
    except Exception as e:
        log_error("Unable to retrieve file contents for {0}:\n{1}".format(item_name, e))

if __name__ == "__main__":
    
    data_file  = open("data.json", "r")
    data = json.load(data_file)
    data_file.close()

    # parameters
    names = ['four_phen', 'musculus', 'neurospora', 'arabidopsis']
    classical_parameters = json.loads(get_item("parameters", "classical-parameters.json").decode('utf8'))

    job_index = int(os.environ['JOB_INDEX'])

    if job_index < 10000:
        gspace_name = names[1]
        initial_node_index = job_index//1000
    elif job_index < 20000:
        gspace_name = names[2]
        initial_node_index = (job_index-10000)//1000
    elif job_index < 30000:
        gspace_name = names[3]
        initial_node_index = (job_index-20000)//1000
    else:
        gspace_name = names[0]
        initial_node_index = (job_index-30000)//1000

    gspace_file = get_item("general-data", data[gspace_name]['g-space']).decode("utf8")
    gspace = nx.parse_gml(gspace_file, label='id')

    phenotypes_file = get_item("general-data", data[gspace_name]["phenotypes"]).decode("utf8")
    phenotypes = phenotypes_file.splitlines()

    genotype_networks =json.loads(get_item("general-data", data[gspace_name]["g-networks"]).decode("utf8"))
    initial_nodes = classical_parameters['initial_nodes'][gspace_name]

    initial_genotype = initial_nodes[initial_node_index]
    max_simulation_time = classical_parameters['max_time'] 
    gamma_c = classical_parameters['transition_rate']
    
    classic_simulation.simulation_cw(gspace, gspace_name, phenotypes, initial_genotype, max_simulation_time, gamma_c)