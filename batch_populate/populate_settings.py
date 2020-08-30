import datetime
settings = {"display_progress": True, 
            "reserve_jobs": True,
            "suppress_errors": True,
            "order": "random"}

def now():
    return datetime.datetime.today().strftime("%Y-%m-%d | %H:%M:%S")