import time

import requests
from cdci_data_analysis.analysis.products import QueryOutput
from cdci_data_analysis.configurer import DataServerConf

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
        excep = Exception()
        
        print('url', self.data_server_url)
        url = self.data_server_url

        for i in range(max_trial):
            try:
                res = requests.get("%s" % (url), params=None)
                print('status_code',res.status_code)
                if res.status_code !=200:
                    no_connection =True
                    raise ConnectionError(f"Backend connection failed: {res.status_code}")
                else:
                    no_connection=False

                    message = 'Connection OK'
                    query_out.set_done(message=message, debug_message='OK')
                    print('-> test connections passed')
                    break
            except Exception as e:
                excep = e
                no_connection = True

            time.sleep(sleep_s)

        if no_connection is True:
            query_out.set_query_exception(excep, 
                                          'no data server connection',
                                          logger=logger)
            raise ConnectionError('Backend connection failed')

        return query_out
    

    def test_has_input_products(self, instrument, logger=None):
        query_out = QueryOutput()
        query_out.set_done('input products check skipped')
        return query_out, []
    
        print('--> test for data availability')

        query_out = QueryOutput()
        streaks = []
        
        url = self.data_server_url + '/api/v1.0/get/checkdata'
        print('url', url)
        
        t1 = instrument.get_par_by_name('T1').value
        t2 = instrument.get_par_by_name('T2').value
        detector = instrument.get_par_by_name('detector').value
        
        res = requests.get("%s" % (url), 
                           params={'t1': t1, 't2': t2, 'detector': detector})
        if res.status_code == 200:
            _o_dict = res.json()
            if _o_dict['output']['ok_flag'] is True:
                #streaks.append(1) # dummy
                query_out.set_done('streak data available')
            else:
                query_out.set_failed('no data available')                
                raise RuntimeError('no data available')
        else:
            excep = RuntimeError('error checking data availability')
            query_out.set_query_exception('error checking data availability', 
                                          message='connection status code: ' + str(res.status_code), 
                                          debug_message=res.text)
            raise excep
            
        return query_out, streaks
    

    def run_query(self,
                  call_back_url = None,
                  run_asynch = True,
                  logger = None,
                  task = None,
                  param_dict = None):
        
        res = None
        message = ''
        debug_message = ''
        query_out = QueryOutput()

        if task is None:
            task=self.task     

        if param_dict is None:
            param_dict=self.param_dict   
        
        if run_asynch:
            param_dict['_async_request_callback'] = call_back_url
            param_dict['_async_request'] = "yes"

        url = '/'.join([self.data_server_url.strip('/'), task.strip('/')])
        res = requests.get(url, params = param_dict)
        if res.status_code == 200:
            if res.json()['data']['exceptions']: #failed nb execution in async 
                except_message = res.json()['data']['exceptions'][0]['ename']+': '+res.json()['data']['exceptions'][0]['evalue']
                query_out.set_failed('Processing failed', 
                                     message=except_message)
                raise RuntimeError(f'Processing failed. {except_message}')
            if task.split('/')[-1] == 'conesearch' and param_dict['do_cone_search']:
                comment = ''
            else:
                comment = 'Please note that RA and Dec parameters are not used.'
            query_out.set_done(message=message, debug_message=str(debug_message),job_status='done',comment=comment)
        elif res.status_code == 201:
            if res.json()['workflow_status'] == 'submitted':
                query_out.set_status(0, message=message, debug_message=str(debug_message),job_status='submitted')
            else:
                query_out.set_status(0, message=message, debug_message=str(debug_message),job_status='progress')
                #this anyway finally sets "submitted", the only status implemented now in "non-integral" dispatcher code
        else:
            try:
                query_out.set_failed('Error in the backend', 
                                 message='connection status code: ' + str(res.status_code), 
                                 extra_message=res.json()['exceptions'][0])
            except:
                query_out.set_failed('Error in the backend', 
                                 message='connection status code: ' + str(res.status_code), 
                                )
            raise RuntimeError('Error in the backend')

        return res, query_out
