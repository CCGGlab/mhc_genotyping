curl 'https://www.internationalgenome.org/api/beta/sample/_search/igsr_samples.tsv' -H 'User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:87.0) Gecko/20100101 Firefox/87.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8' -H 'Accept-Language: en-US,en;q=0.5' --compressed -H 'Content-Type: application/x-www-form-urlencoded' -H 'Origin: https://www.internationalgenome.org' -H 'Connection: keep-alive' -H 'Referer: https://www.internationalgenome.org/data-portal/sample' -H 'Cookie: data_notice_20180523=dismiss' -H 'Upgrade-Insecure-Requests: 1' --data-raw 'json=%7B%22fields%22%3A%5B%22name%22%2C%22sex%22%2C%22biosampleId%22%2C%22populations.code%22%2C%22populations.name%22%2C%22populations.superpopulationCode%22%2C%22populations.superpopulationName%22%2C%22populations.elasticId%22%2C%22dataCollections.title%22%5D%2C%22column_names%22%3A%5B%22Sample+name%22%2C%22Sex%22%2C%22Biosample+ID%22%2C%22Population+code%22%2C%22Population+name%22%2C%22Superpopulation+code%22%2C%22Superpopulation+name%22%2C%22Population+elastic+ID%22%2C%22Data+collections%22%5D%7D'