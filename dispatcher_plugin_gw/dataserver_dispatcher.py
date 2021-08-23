from cdci_data_analysis.analysis.products import QueryOutput
from cdci_data_analysis.configurer import DataServerConf
import requests
import time

class GWDispatcher:
    def __init__(self, instrument=None, param_dict=None, task=None, config=None):
        if config is None:
            config = DataServerConf.from_conf_dict(instrument.data_server_conf_dict)
        self.data_server_url = config.data_server_url
        self.task = task
        self.param_dict = param_dict

    def test_communication(self, max_trial=10, sleep_s=1, logger=None):
        print('--> start test connection')

        query_out = QueryOutput()
        no_connection = True
        debug_message='OK'

        print('url', self.data_server_url)
        url = self.data_server_url

        for i in range(max_trial):
            try:
                res = requests.get("%s" % (url), params=None)
                print('status_code',res.status_code)
                if res.status_code !=200:
                    no_connection =True
                    e = ConnectionError(f"Backend connection failed: {res.status_code}")
                else:
                    no_connection=False

                    message = 'Connection OK'
                    query_out.set_done(message=message, debug_message=str(debug_message))
                    break
            except Exception as e:
                no_connection = True

            time.sleep(sleep_s)

        if no_connection is True:
            message = 'no data server connection'
            debug_message = 'no data server connection'
            connection_status_message = 'no data server connection'

            query_out.set_failed(message,
                                 message='connection_status=%s' % connection_status_message,
                                 logger=logger,
                                 excep=e,
                                 e_message=message,
                                 debug_message=debug_message)

            raise Exception('Connection Error', debug_message)

        print('-> test connections passed')

        return query_out

    def test_has_input_products(self, instrument, logger=None):
        query_out = self.test_communication(logger=logger)
        return query_out, [1] #dummy

    def run_query(self,
                  call_back_url=None,
                  run_asynch = False, #TODO: it should really be True in most cases. To test
                  logger=None,
                  task = None,
                  param_dict=None):
        
        res = None
        message = ''
        debug_message = ''
        query_out = QueryOutput()

        if task is None:
            task=self.task     

        if param_dict is None:
            param_dict=self.param_dict   

        #TODO: handle fail
        url = "%s/%s" % (self.data_server_url, task)
        res = requests.get("%s" % (url), params = param_dict)
        query_out.set_done(message=message, debug_message=str(debug_message),job_status='done')

        return res, query_out