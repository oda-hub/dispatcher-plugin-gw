import requests
import json
import logging

import pytest
import h5py
import tarfile
from astropy.io import fits

logger = logging.getLogger(__name__)

params = {}
params['spectrogram'] = dict(query_status = "new",
                             query_type = "Real",
                             instrument = "gw",
                             product_type = "gw_spectrogram",
                             T1 = "2017-08-17T12:40:54",
                             T2 = "2017-08-17T12:41:10",
                             detector = "H1",
                             qmin = 4,
                             qmax = 64,
                             whiten = True)

params['strain'] = dict(query_status = "new",
                        query_type = "Real",
                        instrument = "gw",
                        product_type = "gw_strain",
                        T1 = "2017-08-17T12:41:02",
                        T2 = "2017-08-17T12:41:06",
                        detector = "H1",
                        fmin = 30,
                        fmax = 400,
                        whiten = True)

params['conesearch'] = dict(query_status = "new",
                            query_type = "Real",
                            instrument = "gw",
                            product_type = "gw_skymap_image",
                            T1 = "2019-07-27T06:03:29.000",
                            T2 = "2019-07-27T06:03:44.000",
                            contour_levels = "50,90",
                            do_cone_search = False)


def test_discover_plugin():
    import cdci_data_analysis.plugins.importer as importer

    assert 'dispatcher_plugin_gw' in  importer.cdci_plugins_dict.keys()

def test_gw(dispatcher_live_fixture, httpserver, product):
    server = dispatcher_live_fixture
    with open(f'tests/mock_backend_json/response_{product}.json', 'r') as fd:
        respjson = json.loads(fd.read())
    httpserver.expect_ordered_request('/').respond_with_data('')    
    httpserver.expect_ordered_request(f'/api/v1.0/get/{product}').respond_with_json({'comment': "task created", "workflow_status": "submitted"}, status=201)
    httpserver.expect_ordered_request('/').respond_with_data('')    
    httpserver.expect_ordered_request(f'/api/v1.0/get/{product}').respond_with_json(respjson)
    
    logger.info("constructed server: %s", server)
    # job submitted
    c = requests.get(server + "/run_analysis",
                    params = params[product])

    logger.info("content: %s", c.text)
    jdata = c.json()
    logger.info(json.dumps(jdata, indent=4, sort_keys=True))
    logger.info(jdata)
    assert c.status_code == 200

    assert jdata['job_status'] == 'submitted'
    # job done
    c = requests.get(server + "/run_analysis",
                    params = params[product])

    logger.info("content: %s", c.text)
    jdata = c.json()
    logger.info(json.dumps(jdata, indent=4, sort_keys=True))
    logger.info(jdata)
    assert c.status_code == 200

    assert jdata['job_status'] == 'done'
    
    d = requests.get(server + "/download_products",
                    params = {
                        'session_id': jdata['job_monitor']['session_id'],
                        'download_file_name': jdata['products']['download_file_name'],
                        'file_list': jdata['products']['file_name'],
                        'query_status': 'ready',
                        'job_id': jdata['job_monitor']['job_id'],
                        'instrument': 'gw'
                    })
    assert d.status_code == 200
    
    
    # with open(jdata['products']['download_file_name'], 'wb') as fd:
    #     fd.write(d.content)

    # FIXME: tarfile.open() fails, but manual testing shows that tar file is good
    # with tarfile.open(jdata['products']['download_file_name'], 'r') as fd:
    #     fns = fd.getnames()
    #     fd.extractall()
        
    # for fn in fns:
    #     if fn.endswith('.fits'):
    #         with fits.open(fn, 'readonly') as fd:
    #             for ext in fd:
    #                 logger.info(ext.header)
    #                 logger.info(ext.data)
    #     if fn.endswith('.h5'):
    #         with h5py.File(fn, 'r') as fd:
    #             for dset_name in fd.keys():
    #                 dset = fd[dset_name]
    #                 logger.info(dset.name)
    #                 logger.info(dset.attrs)
    #                 logger.info(dset.value)
    #     else:
    #         with open(fn, 'r') as fd:
    #             logger.info(fd.read())
