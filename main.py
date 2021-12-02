import classic_simulation
import json
import os
import networkx as nx
import ibm_boto3
from ibm_botocore.client import Config
from ibm_botocore.exceptions import ClientError
from ibmcloudant.cloudant_v1 import CloudantV1, Document
from ibm_cloud_sdk_core.authenticators import IAMAuthenticator
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
    num_simulations = classical_parameters['num_simulations']

    # database connection
    authenticator = IAMAuthenticator(os.environ['CLOUDANT_APIKEY'])

    cloudant = CloudantV1(authenticator=authenticator)
    cloudant.set_service_url(os.environ['CLOUDANT_URL'])

    try:    
        client = cloudant.new_instance()
    except:
        print("Could not create a cloudant client")
    
    for i in range(num_simulations):
        print("Simulation {0}".format(i))
        results = classic_simulation.simulation_cw(gspace, phenotypes, initial_genotype, max_simulation_time, gamma_c)

        simulation: Document = Document()

        simulation.initial_gen_index = results["initial_gen_index"]
        simulation.initial_gen = results["initial_gen"]
        simulation.initial_phen = results["initial_phen"]
        simulation.transition_rate = results["transition_rate"]
        simulation.max_simulation_time = results["max_simulation_time"]
        simulation.total_mutations = results["total_mutations"]
        simulation.computing_time = results["computing_time"]
        simulation.simulation_time = results["simulation_time"]
        simulation.date = results["date"]

        for phen in phenotypes:
            setattr(simulation, 'tau_'+phen, results["tau_"+phen])
            setattr(simulation, 'mutations_'+phen, results["mutations_"+phen])

        try:
            client.post_document(
            db="simulations-cw-"+gspace_name,
            document=simulation
            )
            print("Wrote in Cloudant successfully")

        except:
            print("Unexpected error")
