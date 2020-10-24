import sys
import time
import datajoint as dj

sys.path.append('..')

# Load personal schema and populate settings
from dj_schemas.bvs import BVField, BVScore
from populate_settings import now, settings

dj.config['loglevel'] = 'DEBUG'

def main():
    while True:
        print('############ BVField ############')
        print('{} | BVField.populate()'.format(now()))
        BVField.populate(**settings)

        print('############ BVScore ############')
        print('{} | BVScore.populate()'.format(now()))
        BVScore.populate(**settings)

        time.sleep(5)


if __name__ == "__main__":
    main()