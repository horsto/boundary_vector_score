import sys
import time
import datajoint as dj

sys.path.append('..')

# Load personal schema and populate settings
from dj_schemas.shuffling_bvs import * 
from populate_settings import now, settings

dj.config['loglevel'] = 'DEBUG'

def main():
    while True:
        print('############ ShuffledBVS ############')
        print('{} | ShuffledBVS.populate()'.format(now()))
        ShuffledBVS.populate(**settings)

        time.sleep(15)


if __name__ == "__main__":
    main()