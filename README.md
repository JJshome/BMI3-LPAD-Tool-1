# LAMP Primer Auto Design (LPAD) Tool

## 개요
LAMP Primer Auto Design (LPAD) Tool은 Loop-mediated Isothermal Amplification (LAMP)을 위한 프라이머 설계를 자동화하는 도구입니다. 이 도구는 녹는 온도(Tm), GC 함량, 터미널 ΔG, 프라이머 거리, 2차 구조와 같은 주요 특성에 기반하여 프라이머에 점수를 매깁니다. 프라이머 특이성은 사용자 정의 참조 게놈을 사용한 BLAST 정렬을 통해 평가됩니다. 그런 다음 Random Forest 모델을 사용하여 프라이머 품질을 예측하고, 프라이머 세트의 순위를 매깁니다. 실제 및 시뮬레이션 데이터 세트에 대한 평가는 높은 모델 정확도를 보여줍니다. LPAD는 기존 도구에 비해 더 효율적이고 포괄적인 LAMP 프라이머 설계 솔루션을 제공합니다.

## 목차
- [환경 설정](#환경-설정)
- [전체 파이프라인](#전체-파이프라인)
- [구체적인 구현](#구체적인-구현)
- [사용법](#사용법)
- [파일 구조](#파일-구조)

## 환경 설정
환경 설정은 `environment.yml` 파일에 제공됩니다.

### 1. 운영 체제
이 도구를 실행하기 위해 Linux 또는 MacOS 사용을 강력히 권장합니다.

### 2. 의존성
- python=3.8
- primer3-py=2.0.3
- biopython=1.78
- pandas=2.0.3
- blast=2.16.0
- joblib=1.4.2
- scipy=1.10.1
- scikit-learn=1.3.0

## 전체 파이프라인

### 1. 설계 목적：
LAMP 증폭에는 forward inner primer(FIP), backward inner primer(BIP), 두 개의 outer primer(F3 및 B3), 그리고 두 개의 loop primer(LF 및 LB)와 같은 6개의 프라이머가 필요합니다. 이 프라이머들은 특정 서열 영역을 표적으로 하는 데 사용됩니다. 그 중 루프 프라이머는 선택적입니다. 이 프로젝트는 주로 forward inner primer(FIP), backward inner primer(BIP), 두 개의 outer primer(F3 및 B3)를 자동으로 설계하고, 고품질 F1, F2, F3; B1, B2, B3 프라이머를 자동으로 출력하는 데 중점을 둡니다.

### 2. 도구 워크플로우：
#### 저장소 클론
```
git clone https://github.com/JJshome/BMI3-LPAD-Tool-1
cd BMI3-LPAD-Tool-1
```

#### 환경 설정:
컴퓨터에 Conda가 설치되어 있는지 확인하세요. 필요한 의존성을 설치하려면 제공된 environment.yml 파일을 Conda와 함께 사용하세요:
```
conda env create -f environment.yml
```
이렇게 하면 프로젝트에 필요한 모든 의존성을 갖춘 새로운 Conda 환경이 생성됩니다. 이 프로그램을 사용하려면 현재 프로그램의 홈 디렉토리로 이동하고 conda 환경을 활성화하세요:
```
conda activate lpad
```

#### 입력 데이터：
도구 입력에는 타겟 시퀀스(증폭하려는 영역)의 정렬 정보와 배경 게놈 정보(증폭하지 않으려는 영역)가 포함됩니다.
```
# 기본 경로 사용 (*** 참고!!! 이 방법은 배경 게놈의 fasta 파일 'hg38.fa'가 './data/resource' 폴더에 미리 배치되어 있어야 합니다 ***)
python LPAD.py
# BMI3-LPAD-Tool 폴더 아래에 사용자 지정 경로를 지정할 수도 있습니다
python LPAD.py -i /path/to/your/input.fasta -r /path/to/your/reference.fasta -o /path/to/output/directory
```
개별 LAMP 프라이머에 대한 매개변수를 지정하려면 './configs' 폴더의 해당 프라이머 파일에 대한 매개변수를 수정하세요.

#### 시퀀스 정렬 및 프라이머 생성:
프라이머는 특정 제약(길이, 간격)에 따라 대상 시퀀스에서 생성되며, 배경(오염) 게놈에 결합하는 것을 방지하기 위해 특이성이 엄격하게 검사됩니다.

#### 제약 조건 정의 및 점수 기준:
각 후보 프라이머는 제약 조건(Tm, GC 함량, 자유 에너지 등)에 대해 평가되며, 이러한 조건을 얼마나 잘 충족하는지에 따라 점수가 할당됩니다.
2차 구조 예측 및 이합체 형성 필터링이 추가되어 있어, 점수 정확도를 더욱 향상시키고 잠재적으로 문제가 될 수 있는 프라이머를 제거합니다.

#### 결과 출력:
이 도구는 사용자 선택을 위해 가장 높은 점수에서 가장 낮은 점수 순으로 정렬된 프라이머 조합 목록을 출력합니다. 출력 파일은 './data/output/Final_score/'에 저장되며, 다른 중간 파일들은 './data/output/Intermediate_file/'에 저장됩니다.

## 구체적인 구현:
### 알고리즘의 실현 가능성:
#### TM:
- Tm은 Nearest-Neighbor 방법을 사용하여 추정됩니다. 이 방법은 현재 실제 값에 가장 가까운 근사값을 제공하는 방법으로 간주됩니다.
- 계산된 Tm은 염 농도 및 올리고 농도와 같은 실험 조건의 영향을 받으므로, Tm은 고정된 실험 조건(올리고 농도 0.1 µM, 나트륨 이온 농도 50 mM, 마그네슘 이온 농도 4 mM)에서 계산되는 것이 선호됩니다.
- 각 영역의 Tm은 F1c와 B1c의 경우 약 65°C(64 - 66°C), F2, B2, F3, B3의 경우 약 60°C(59 - 61°C), 루프 프라이머의 경우 약 65°C(64 - 66°C)로 설계됩니다.

#### GC%:
- 프라이머는 GC 함량이 약 40%에서 65% 사이가 되도록 설계됩니다.
- GC 함량이 50%에서 60% 사이인 프라이머는 상대적으로 좋은 프라이머를 제공하는 경향이 있습니다.

#### ΔG:
- 프라이머의 끝은 DNA 합성의 시작점 역할을 하므로 일정 수준의 안정성이 있어야 합니다. F2/B2, F3/B3, LF/LB의 3' 말단과 F1c/B1c의 5' 말단은 자유 에너지가 –4 kcal/mol 이하가 되도록 설계됩니다. 증폭 후 F1c의 5' 말단은 F1의 3' 말단에 해당하므로 안정성이 중요합니다.

#### 프라이머 거리:
- (Primer Explorer 문서 참조)

#### 2차 구조 예측 및 이합체 형성:
- 서열 내 상보적 서열(hairpin)이 없고 서열 간 상보적 서열(dimer)이 없는지 고려합니다.
- Hairpin 설계: 1. Hairpin 구조의 상보적 세그먼트 길이는 6-12 bp 사이여야 합니다. 2. 루프 세그먼트의 길이는 4-8 bp 범위 내여야 합니다.
- Dimer 설계: 서열 간에 8-16bp 상보적 영역이 있는지 확인합니다.

## 사용법

### 기본 사용법
```bash
python LPAD.py
```

### 고급 사용법
```bash
python LPAD.py -i input.fasta -r reference.fasta -o output_directory
```

### 배치 모드
```bash
python LPAD.py -b -i input_directory -r reference.fasta -o output_base_directory
```

### 명령줄 옵션
- `-i, --input_file` : 입력 FASTA 파일 또는 배치 모드의 디렉토리 경로
- `-r, --ref_file` : 참조 FASTA 파일 경로
- `-o, --output_dir` : 출력이 저장될 디렉토리
- `-b, --batch` : 입력 파일 디렉토리를 사용한 배치 모드로 실행
- `-p, --pattern` : 배치 모드용 파일 패턴 (기본값: *.fasta)
- `-c, --clean` : 실행 전 출력 디렉토리 정리
- `-v, --verbose` : 자세한 출력 인쇄
- `--version` : 버전 정보 표시

## 파일 구조

```
BMI3-LPAD-Tool-1/
├── LPAD.py               # 메인 실행 스크립트
├── README.md             # 프로젝트 문서
├── environment.yml       # Conda 환경 설정
├── batch_process.py      # 배치 처리 모듈
├── calc_features_multi.py  # 다중 프라이머 특성 계산
├── calc_features_single.py # 단일 프라이머 특성 계산
├── final_output.py         # 최종 결과 처리
├── heuristic_primer_design.py  # 프라이머 설계 알고리즘
├── primer_score.py       # 프라이머 점수 계산
├── specificity.py        # 특이성 검사
├── configs/              # 구성 파일
│   ├── inner_primer.json   # 내부 프라이머 매개변수
│   ├── loop_primer.json    # 루프 프라이머 매개변수
│   ├── middle_primer.json  # 중간 프라이머 매개변수
│   └── outer_primer.json   # 외부 프라이머 매개변수
├── core/                 # 코어 알고리즘
│   ├── __init__.py
│   ├── combination.py    # 프라이머 조합
│   ├── mutant.py         # 변이 처리
│   └── overlap.py        # 중복 처리
├── data/                 # 데이터 디렉토리
│   ├── input/            # 입력 파일
│   │   └── example.fasta  # 예제 입력
│   ├── output/           # 출력 디렉토리
│   │   ├── Final_score/    # 최종 결과
│   │   └── Intermediate_file/  # 중간 파일
│   ├── resource/         # 리소스 파일
│   │   └── hg38.fa        # 참조 게놈 (매우 큰 파일)
│   └── rfmodel_data/     # 모델 훈련 데이터
├── model/                # 모델 디렉토리
│   └── random_forest_model.pkl  # 훈련된 모델
└── utils/                # 유틸리티 모듈
    ├── __init__.py
    ├── config_loader.py  # 구성 로더
    ├── logging_utils.py  # 로깅 유틸리티
    ├── oligos.py         # 올리고 처리
    ├── penalty.py        # 페널티 계산
    ├── pipeline.py       # 파이프라인 실행
    ├── runprimer3.py     # Primer3 실행
    ├── sequence_utils.py # 시퀀스 유틸리티
    └── validation.py     # 입력 검증
```
