#! /usr/bin/env python

import sys, re
import subprocess

class Storage(object):

    def __init__(self, region = None, zone = None, dryrun = False):

        """
        if region is None and zone is None:
            print ("region or zone should not be None.")
            sys.exit(1)
        """
        self.provider = None
        self.region = region
        self.zone = zone
        self.dryrun = dryrun

    def upload(self, local_file_path, storage_path, create_bucket = False):
        print ("=== Storage upload ===")
        print ("local_file_path: %s" % (local_file_path))
        print ("storage_path: %s" % (storage_path))
        print ("create_bucket: %s" % (create_bucket))
        
        if self.dryrun:
            return
            
        if storage_path.startswith("s3://"):
            self.__upload_to_aws(local_file_path, storage_path, create_bucket)
        elif storage_path.startswith("gs://"):
            self.__upload_to_gcs(local_file_path, storage_path, create_bucket)
        else:
            raise ValueError ("Storage path should starts with 's3://' or 'gs://'.")


    def __get_bucket_name_from_storage_path(self, storage_path):

        if storage_path.startswith("s3://"):
            return re.sub(r'^s3://', '', storage_path).split('/')[0]
        elif storage_path.startswith("gs://"):
            return re.sub(r'^gs://', '', storage_path).split('/')[0]
        else:
            raise ValueError ("Storage path should starts with 's3://' or 'gs://'.")


    def __print_error (self, message):
        if sys.version_info.major == 2:
            print >> sys.stderr, message    
        else:
            #print (message, file = sys.stderr)
            print (message)
            
    # #####################
    # upload
    # #####################
    def __upload_to_aws(self, local_file_path, storage_path, create_bucket = False):

        # Extract AWS S3 bucket name of the storage_path argument
        target_bucket_name = self.__get_bucket_name_from_storage_path(storage_path)

        if not self.__check_bucket_exist_aws(target_bucket_name):
            if create_bucket:
                self.__print_error ("bucket %s is not exist." % (target_bucket_name))
                self.__create_bucket_aws(bucket_name = target_bucket_name)
            else:
                self.__print_error ("No bucket: " + target_bucket_name)
                sys.exit(1)

        # Upload to the S3
        import boto3
        storage_file_name = '/'.join(re.sub(r'^s3://', '', storage_path).split('/')[1:])
        boto3.client("s3").upload_file(local_file_path, target_bucket_name, storage_file_name) 


    def __upload_to_gcs(self, local_file_path, storage_path, create_bucket = False):

        # Extract Google Cloud Storage bucket name if the storage_path argument
        target_bucket_name = self.__get_bucket_name_from_storage_path(storage_path)

        if not self.__check_bucket_exist_gcs(target_bucket_name):

            if create_bucket:
                proc = subprocess.Popen(["gsutil", "mb", 'gs://' + target_bucket_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
                out, err = proc.communicate()
                if proc.returncode == 1:
                    if "ServiceException" in err:
                        self.__print_error ("Bucket " + target_bucket_name + " already exists.")
                        sys.exit(1)
                    else:
                        self.__print_error ("An Error happend while creating Buckt " + target_bucket_name + ".")
                        sys.exit(1)
            else:
                self.__print_error ("No bucket: " + target_bucket_name)
                sys.exit(1)

        # Upload to the Google Cloud Storage
        proc = subprocess.Popen(["gsutil", "cp", local_file_path, storage_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        proc.communicate()


    # #####################
    # check bucket exists
    # #####################
    def __check_bucket_exist_aws(self, bucket_name):

        import boto3

        s3 = boto3.client("s3")
        response = s3.list_buckets()
        available_bucket_names = [bucket["Name"] for bucket in response["Buckets"]]

        if bucket_name in available_bucket_names:
            return True
        else:
            return False


    def __check_bucket_exist_gcs(self, bucket_name):

        proc = subprocess.Popen(["gsutil", "ls"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        available_bucket_names = [re.sub(r'^gs://', '', x).rstrip('/') for x in proc.communicate()[0].split('\n')]

        if bucket_name in available_bucket_names:
            return True
        else:
            return False


    # #####################
    # create bucket
    # #####################
    def __create_bucket_aws(self, bucket_name):
        
        proc = subprocess.Popen(["aws", "s3", "mb", "s3://" + bucket_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = proc.communicate()
        if proc.returncode != 0:
            self.__print_error (err)
            sys.exit(1)
                
    def __create_bucket_gcs(self, bucket_name):

        proc = subprocess.Popen(["gsutil", "mb", 'gs://' + bucket_name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = proc.communicate()
        if proc.returncode == 1:
            if "ServiceException" in err:
                self.__print_error ("Bucket " + bucket_name + " already exists.")
                sys.exit(1)
            else:
                self.__print_error ("An Error happend while creating Buckt " + bucket_name + ".")
                sys.exit(1)


if __name__ == "__main__":

    storage = Storage()
    storage.upload("test.txt", "gs://friend1_upload_test/test.txt", True)  

