import base64
import hashlib
import json
import os
import re
import uuid
from copy import deepcopy
from time import sleep

import requests


# Although this was originally intended to be async, its actual usage is fully
# synchronous: the server either retrieves the SMARTS image from the cache or
# generates it on demand, returning the result immediately and updating the cache.
# Since this will run inside a standard for-loop over a dataframe in a notebook,
# a synchronous implementation is more appropriate.


class SmartsRenderWorker:
    def __init__(self):
        self.endpoint = "https://api.smarts.plus/smartsView/"
        self.headers = {
            "Content-Type": "application/json",
            "X-API-Key": os.environ["SMARTS_API_KEY"],
        }
        self._body_template = {
            "job_id": "your_job_id",
            "query": {
                "smarts": "SMARTS_pattern",
                "parameters": {"file_format": "svg"},
            },
        }

    def generate_smarts(self, smarts: str, format: str):
        job_id, r = self.submit_job(smarts, format)

        # If the server returns a 200 from submission,
        # the results are already gone on the server
        # *** requests.exceptions.HTTPError: 404 Client Error:
        # Not Found for url: https://api.smarts.plus/smartsView/?job_id=myjob1111
        # Hence we return them as we get them
        if r.status_code == 200:
            return SmartsRenderWorker.process_output(r)

        # If the response is not 200, however,
        # we need to do polling until they're ready
        return self.fetch_results(job_id)

    def process_output(response: requests.models.Response):
        r_json = response.json()
        if r_json.get("result").get("image"):
            return r_json["result"]["image"]
        elif r_json.get("query").get("result"):
            if "error" in r_json["query"]["result"].keys():
                return None

    def submit_job(self, smarts: str, format: str):
        job_id = str(uuid.uuid1())
        payload = self._prepare_payload(job_id, smarts, format)
        r = requests.post(
            self.endpoint,
            headers=self.headers,
            data=json.dumps(payload),
        )
        r.raise_for_status()
        return job_id, r

    def _prepare_payload(self, job_id: str, smarts: str, format: str):
        payload = deepcopy(self._body_template)
        payload["job_id"] = job_id
        payload["query"]["smarts"] = smarts
        payload["query"]["parameters"]["file_format"] = format
        return payload

    def fetch_results(self, job_id):
        while True:
            r = requests.get(f"{self.endpoint}?job_id={job_id}")
            r.raise_for_status()
            status = r.status_code
            if status == 200:
                return SmartsRenderWorker.process_output(r)
            elif status == 202:
                sleep(1)
                continue

    async def close(self):
        await self.client.aclose()


class SmartsRenderServer:
    def __init__(self):
        self.cache_dir = os.path.join(os.path.dirname(__file__), "cache")
        self._error_svg = open(
            os.path.join(os.path.dirname(__file__), "assets", "error.svg")
        ).read()
        self.worker = SmartsRenderWorker()

    def generate_smarts_base64_svg(self, smarts: str, width=100, height=100):
        full_svg = self.generate_smarts_image(smarts, format="svg")
        contents = re.sub(r"<\?xml.*?\?>", "", full_svg)
        b64 = base64.b64encode(contents.encode("utf-8")).decode("utf-8")
        return f'<img src="data:image/svg+xml;base64,{b64}" width="{width}" height="{height}">'  # noqa

    def generate_smarts_image(self, smarts: str, format: str = "svg"):
        filename = SmartsRenderServer._make_hash(smarts) + f".{format}"
        path = os.path.join(self.cache_dir, filename)
        image_str = SmartsRenderServer.load_if_exists(path)
        if not image_str:
            image_str = self.worker.generate_smarts(smarts, format=format)
            if not image_str:
                image_str = self._error_svg
            SmartsRenderServer.write_to_cache(image_str, path)
        return image_str

    def _make_hash(smarts: str):
        """Deterministic 16-char short hash."""
        digest = hashlib.md5(smarts.encode()).digest()
        return base64.urlsafe_b64encode(digest)[:16].decode()

    def load_if_exists(path):
        return open(path).read() if os.path.isfile(path) else None

    def write_to_cache(image_str: str, path):
        with open(path, "w") as w:
            w.write(image_str)


if __name__ == "__main__":
    svg_smarts = SmartsRenderServer().generate_smarts("[SX4](=O)(=O)[OX2H]")
