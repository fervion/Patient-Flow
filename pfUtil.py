from random import random

"""Utility Functions used by patient flow"""
#draws from a weighted distribution where distribution is a dictionary with choices to weights.  weights are not normalized
def draw_from_distribution(distribution):
    sumWeights = sum(distribution.values())
    randNum = random()*sumWeights
    culNum = 0
    for key, value in distribution.iteritems():
        culNum += value
        if randNum <= culNum:
            return key
        
def isWeekday(date):
    weekends = (5,6)
    if date.weekday() in weekends:
        return False
    else:
        return True
    
class InvalidLayoutError(Exception):
    pass
