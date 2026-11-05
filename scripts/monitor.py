#!/usr/bin/env python3
"""Weekly PubMed discovery; metadata only, no eligibility inference."""
import json, os, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; PROJECTS=DATA/'projects'; PROJECTS.mkdir(parents=True,exist_ok=True)
GLP='("Glucagon-Like Peptide-1 Receptor Agonists"[Mesh] OR semaglutide[tiab] OR liraglutide[tiab] OR tirzepatide[tiab] OR dulaglutide[tiab] OR exenatide[tiab] OR retatrutide[tiab] OR "GLP-1"[tiab] OR incretin*[tiab])'
EX='(exercise[tiab] OR training[tiab] OR physical activity[tiab] OR resistance training[tiab])'
TERMS={
1:f'{GLP} AND {EX} AND (muscle[tiab] OR strength[tiab] OR physical function[tiab])',
2:f'{GLP} AND (VO2[tiab] OR peak oxygen[tiab] OR cardiorespiratory fitness[tiab] OR exercise capacity[tiab])',
3:'(exercise[tiab] OR training[tiab]) AND (adipose tissue[tiab]) AND (inflammation[tiab] OR vascularization[tiab] OR fibrosis[tiab])',
4:f'{GLP} AND (body composition[tiab] OR lean mass[tiab] OR fat mass[tiab]) AND (longitudinal[tiab] OR follow-up[tiab] OR time course[tiab])',
5:f'{GLP} AND (sarcopenia[tiab] OR grip strength[tiab] OR muscle function[tiab] OR physical performance[tiab])',
6:f'{GLP} AND (ventricular function[tiab] OR hemodynamic*[tiab] OR exercise tolerance[tiab] OR echocardiograph*[tiab])',
7:f'{GLP} AND (exercise[tiab] OR protein[tiab] OR nutrition[tiab]) AND (muscle[tiab] OR lean mass[tiab])',
8:f'{GLP} AND (cognition[tiab] OR cognitive function[tiab] OR memory[tiab])',
9:f'{GLP} AND (Alzheimer*[tiab] OR mild cognitive impairment[tiab])',
10:f'{GLP} AND (dementia[tiab] OR neurodegenerative disease[tiab])',
11:f'{GLP} AND (cancer[tiab] OR neoplasm*[tiab] OR malignan*[tiab])',
12:f'{GLP} AND (thyroid cancer[tiab] OR medullary thyroid[tiab] OR thyroid neoplasm*[tiab])',
13:f'{GLP} AND (alcohol[tiab] OR drinking[tiab] OR alcohol use disorder[tiab])',
14:f'{GLP} AND (smoking[tiab] OR tobacco[tiab] OR nicotine[tiab])',
15:f'{GLP} AND (opioid*[tiab] OR stimulant*[tiab] OR cocaine[tiab] OR methamphetamine[tiab])',
16:f'{GLP} AND (autonomic[tiab] OR orthostatic[tiab] OR heart rate variability[tiab] OR baroreflex[tiab])',
17:f'{GLP} AND (sleep apnea[tiab] OR sleep apnoea[tiab] OR obstructive sleep[tiab])',
18:'(tirzepatide[tiab] OR retatrutide[tiab] OR dual agonist[tiab] OR triple agonist[tiab]) AND (obesity[tiab] OR diabetes[tiab])',
19:f'{GLP} AND (discontinuation[tiab] OR withdrawal[tiab] OR weight regain[tiab] OR cessation[tiab])',
20:f'{GLP} AND (predictor*[tiab] OR responder*[tiab] OR heterogeneity[tiab] OR treatment response[tiab])',
21:'(adipose tissue[tiab] OR adipose fibrosis[tiab]) AND (metabolic dysfunction[tiab] OR insulin resistance[tiab]) AND (histolog*[tiab] OR biops*[tiab])',
22:'(brown adipose tissue[tiab] OR brown fat[tiab] OR thermogenesis[tiab]) AND (metabolic health[tiab] OR cardiometabolic[tiab] OR obesity[tiab])',
23:f'{GLP} AND (kidney[tiab] OR renal[tiab] OR albuminuria[tiab] OR eGFR[tiab])',
24:f'{GLP} AND (MASLD[tiab] OR MASH[tiab] OR liver fibrosis[tiab] OR hepatic steatosis[tiab])',
25:f'{GLP} AND (access[tiab] OR equity[tiab] OR affordability[tiab] OR availability[tiab])',
26:f'{GLP} AND (cost effectiveness[tiab] OR cost-effectiveness[tiab] OR economic evaluation[tiab])',
27:f'{GLP} AND (adherence[tiab] OR persistence[tiab] OR discontinuation[tiab])',
28:f'{GLP} AND (quality of life[tiab] OR patient reported outcome*[tiab] OR wellbeing[tiab])'
}
TITLES=['Exercício + GLP-1','GLP-1 e capacidade cardiorrespiratória','Exercício e tecido adiposo','Composição corporal com incretinas','GLP-1 e sarcopenia','GLP-1 e função cardiovascular','Intervenções para preservar massa muscular','GLP-1 e cognição','GLP-1 e doença de Alzheimer','GLP-1 e risco de demência','GLP-1 e câncer','GLP-1 e câncer de tireoide','GLP-1 e álcool','GLP-1 e tabagismo','GLP-1 e dependência de opioides/estimulantes','GLP-1 e sistema nervoso autônomo','GLP-1 e apneia obstrutiva do sono','Novos agonistas múltiplos','Interrupção de GLP-1','Resposta heterogênea ao GLP-1','Tecido adiposo e saúde metabólica','Tecido adiposo termogênico','GLP-1 e função renal','GLP-1 e doença hepática metabólica','GLP-1 e desigualdade de acesso','Custo-efetividade dos GLP-1','GLP-1 e adesão','GLP-1 e qualidade de vida']
EMAIL=os.getenv('NCBI_EMAIL',''); TOOL='glp1_research_observatory'; LIMIT=int(os.getenv('PUBMED_MAX_PER_TOPIC','100'))
def request(endpoint,params):
    params={**params,'tool':TOOL}
    if EMAIL:params['email']=EMAIL
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'?'+urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':TOOL+'/1.0'}),timeout=40) as r:return r.read()
        except Exception:
            if attempt==3:raise
            time.sleep(2**attempt)
def val(e,path):
    node=e.find(path)
    return ' '.join(''.join(node.itertext()).split()) if node is not None else ''
def ris(record):
    fields=['TY  - JOUR','ID  - '+record['pmid'],'TI  - '+record['title']]
    for author in record['authors']:fields.append('AU  - '+author)
    if record['journal']:fields.append('JO  - '+record['journal'])
    if record['year']:fields.append('PY  - '+record['year'])
    if record['doi']:fields.append('DO  - '+record['doi'])
    if record['abstract']:fields.append('AB  - '+record['abstract'].replace('\n',' '))
    fields.extend(['UR  - '+record['url'],'ER  - ',''])
    return '\n'.join(fields)+'\n'
def parse_article(a):
    pmid=val(a,'./MedlineCitation/PMID'); art=a.find('./MedlineCitation/Article')
    if not pmid or art is None:return None
    doi=''
    for idn in a.findall('.//ArticleId'):
        if idn.attrib.get('IdType')=='doi':doi=(idn.text or '').strip();break
    authors=[]
    for au in art.findall('./AuthorList/Author'):
        name=' '.join(x for x in [val(au,'LastName'),val(au,'ForeName')] if x)
        if name:authors.append(name)
    abstract=' '.join(''.join(x.itertext()).strip() for x in art.findall('./Abstract/AbstractText'))
    year=val(art,'./Journal/JournalIssue/PubDate/Year') or val(art,'./Journal/JournalIssue/PubDate/MedlineDate')[:4]
    pubtypes=[(x.text or '').strip() for x in art.findall('./PublicationTypeList/PublicationType')]
    return {'pmid':pmid,'title':val(art,'ArticleTitle'),'authors':authors,'journal':val(art,'./Journal/Title'),'year':year,'doi':doi,'abstract':abstract,'publication_types':pubtypes,'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/','screening_status':'não avaliado'}
summary=[]; errors=[]
for idx in range(1,29):
    folder=PROJECTS/f'{idx:02d}';folder.mkdir(exist_ok=True)
    query=TERMS[idx]+' AND ("2025/01/01"[Date - Publication] : "3000"[Date - Publication])'
    try:
        search=ET.fromstring(request('esearch.fcgi',{'db':'pubmed','term':query,'retmode':'xml','retmax':str(LIMIT),'sort':'pub date'}))
        ids=[x.text for x in search.findall('.//IdList/Id') if x.text]; total=int(val(search,'Count') or 0)
        oldpath=folder/'records.json';old=json.loads(oldpath.read_text()) if oldpath.exists() else []
        byid={r['pmid']:r for r in old}
        if ids:
            for off in range(0,len(ids),50):
                result=ET.fromstring(request('efetch.fcgi',{'db':'pubmed','id':','.join(ids[off:off+50]),'rettype':'xml','retmode':'xml'}))
                for a in result.findall('./PubmedArticle'):
                    rec=parse_article(a)
                    if rec:
                        if rec['pmid'] in byid:rec['screening_status']=byid[rec['pmid']].get('screening_status','não avaliado')
                        byid[rec['pmid']]=rec
                time.sleep(.38)
        records=sorted(byid.values(),key=lambda r:(r.get('year',''),int(r['pmid'])),reverse=True)
        oldpath.write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        (folder/'candidates.ris').write_text(''.join(ris(r) for r in records),encoding='utf-8')
        (folder/'search.json').write_text(json.dumps({'project_id':idx,'title':TITLES[idx-1],'query':query,'last_run_utc':datetime.now(timezone.utc).isoformat(),'pubmed_total_matching':total,'retrieved_in_latest_run':len(ids),'stored_unique_records':len(records),'retrieval_limit_per_run':LIMIT,'note':'Candidatos automáticos não equivalem a artigos incluídos. Estratégia exploratória, requer validação por revisores.'},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        summary.append({'id':idx,'title':TITLES[idx-1],'total':total,'stored':len(records),'new':len(byid)-len(old),'folder':f'data/projects/{idx:02d}/'})
        print(f'{idx:02d}: {len(records)} stored / {total} matching')
    except Exception as exc:
        errors.append({'id':idx,'error':str(exc)[:300]});print(f'ERROR {idx}: {exc}',flush=True)
    time.sleep(.4)
(DATA/'latest.json').write_text(json.dumps({'last_run_utc':datetime.now(timezone.utc).isoformat(),'source':'PubMed (NCBI E-utilities)','summary':summary,'errors':errors,'limitations':'Triagem automatizada apenas por correspondência de busca; não determina elegibilidade, não recupera PDFs nem busca Embase/CENTRAL.'},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
if errors:raise SystemExit(f'{len(errors)} projects failed; successful projects were saved')
