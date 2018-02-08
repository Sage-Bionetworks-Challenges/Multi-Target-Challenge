# Use rpy2 if you have R scoring functions
# import rpy2.robjects as robjects
# import os
# filePath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'getROC.R')
# robjects.r("source('%s')" % filePath)
# AUC_pAUC = robjects.r('GetScores')
##-----------------------------------------------------------------------------
##
## challenge specific code and configuration
##
##-----------------------------------------------------------------------------
import os
import yaml
## A Synapse project will hold the assetts for your challenge. Put its
## synapse ID here, for example
## CHALLENGE_SYN_ID = "syn1234567"
CHALLENGE_SYN_ID = "syn8404040"

## Name of your challenge, defaults to the name of the challenge's project
CHALLENGE_NAME = "Multi-kinase drug inhibitor prediction DREAM challenge"

## Synapse user IDs of the challenge admins who will be notified by email
## about errors in the scoring script
ADMIN_USER_IDS = [3334658]


## Each question in your challenge should have an evaluation queue through
## which participants can submit their predictions or models. The queues
## should specify the challenge project as their content source. Queues
## can be created like so:
##   evaluation = syn.store(Evaluation(
##     name="My Challenge Q1",
##     description="Predict all the things!",
##     contentSource="syn1234567"))
## ...and found like this:
##   evaluations = list(syn.getEvaluationByContentSource('syn3375314'))
## Configuring them here as a list will save a round-trip to the server
## every time the script starts and you can link the challenge queues to
## the correct scoring/validation functions.  Predictions will be validated and 

def validate_func(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Validate submission
    ## MUST USE ASSERTION ERRORS!!! 
    
    # Check file formart
    assert submission.filePath.endswith(".yml") or submission.filePath.endswith(".yaml"),"Submission file must be in yaml file format."
    # Validate/Load file
    with open(submission.filePath,'r') as f:
       try:
           data_loaded = yaml.load(f)
       except yaml.YAMLError:
            assert False,"Your file cannot be parsed. Please check the syntax."
    # Check if required fields are filled
    ## Prediction Methods
    preditction_method = data_loaded['Prediction Methods']
    assert preditction_method and preditction_method != "Replace this text with your prediction methods text, written as for submission to Nature Chemical Biology.",'"Prediction Methods" is required.'
    ## Rationale
    rationale = data_loaded['Rationale']
    rationale_p1 = rationale['Why is your approach innovative?']
    rationale_p2 = rationale['Why will your approach be generalizable?']
    assert rationale_p1 and rationale_p1 != "Replace this text with explanation of why your approach is innovative.",'"Why is your approach innovative?" requires an answer.'
    assert rationale_p2 and rationale_p2 != "Replace this text with explanation of why your approach is generalizable.",'"Why will your approach be generalizable?" requires an answer.'
    ## Solutions
    solution_template = {"ZINC ID": "Your ZINC ID here",
                         "VENDOR ID": "The ID provided for this compound from the vendor",
                         "SMILES string": "Your SMILES string here",
                         "VENDOR NAME": "The name of a vendor from which the compound can be purchased",
                         "Explanation of chemical novelty": "Your explanation here."}
    ## Problem 1 - Solution 1
    p1s1 = data_loaded['Problem 1'][0]['Solution 1']
    for question, answer in p1s1.iteritems():
            assert answer and solution_template[question] != answer,"%s is required for Problem 1 - Solution 1." % question
    ## Problem 1 - Solution 1
    p2s1 = data_loaded['Problem 2'][0]['Solution 1']
    for question, answer in p2s1.iteritems():
            assert answer and solution_template[question] != answer,"%s is required for Problem 2 - Solution 1." % question
    ## Only assertion errors will be returned to participants, all other errors will be returned to the admin
    return(True,"Passed Validation")

def score1(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(5)

def score2(submission, goldstandard_path):
    ##Read in submission (submission.filePath)
    ##Score against goldstandard
    return(score1, score2, score3)

evaluation_queues = [
    {
        'id':9606273,
        #'scoring_func':score1,
        'validation_func':validate_func,
        'goldstandard_path':'path/to/sc1gold.txt'
    }
]
evaluation_queue_by_id = {q['id']:q for q in evaluation_queues}


## define the default set of columns that will make up the leaderboard
LEADERBOARD_COLUMNS = [
    dict(name='objectId',      display_name='ID',      columnType='STRING', maximumSize=20),
    dict(name='userId',        display_name='User',    columnType='STRING', maximumSize=20, renderer='userid'),
    dict(name='entityId',      display_name='Entity',  columnType='STRING', maximumSize=20, renderer='synapseid'),
    dict(name='versionNumber', display_name='Version', columnType='INTEGER'),
    dict(name='name',          display_name='Name',    columnType='STRING', maximumSize=240),
    dict(name='team',          display_name='Team',    columnType='STRING', maximumSize=240)]

## Here we're adding columns for the output of our scoring functions, score,
## rmse and auc to the basic leaderboard information. In general, different
## questions would typically have different scoring metrics.
leaderboard_columns = {}
for q in evaluation_queues:
    leaderboard_columns[q['id']] = LEADERBOARD_COLUMNS + [
        dict(name='score',         display_name='Score',   columnType='DOUBLE'),
        dict(name='rmse',          display_name='RMSE',    columnType='DOUBLE'),
        dict(name='auc',           display_name='AUC',     columnType='DOUBLE')]

## map each evaluation queues to the synapse ID of a table object
## where the table holds a leaderboard for that question
leaderboard_tables = {}


def validate_submission(evaluation, submission):
    """
    Find the right validation function and validate the submission.

    :returns: (True, message) if validated, (False, message) if
              validation fails or throws exception
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    validated, validation_message = config['validation_func'](submission, config['goldstandard_path'])

    return True, validation_message

'''
def score_submission(evaluation, submission):
    """
    Find the right scoring function and score the submission

    :returns: (score, message) where score is a dict of stats and message
             is text for display to user
    """
    config = evaluation_queue_by_id[int(evaluation.id)]
    score = config['scoring_func'](submission, config['goldstandard_path'])
    #Make sure to round results to 3 or 4 digits
    return (dict(score=score), "You did fine!")
'''

