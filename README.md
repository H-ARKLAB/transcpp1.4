# transcpp_1.4

DNA에 결합한 두 TF간의 거리를 3차원 공간 상에서 계산한 새로운 distance function을 포함한 transcpp https://github.com/kennethabarr/transcpp

# 설치법
INSTALL.TXT 참조

transcpp는 xtransc의 업그레이드 버전으로 속도가 많이 향상된 버전입니다. 
2016년 2월 기준 xml파일은 형식이 달라졌으며, unfold 또한 기능이 많이 달라졌습니다. 

프로그램 설치는 conda를 이용해도 무방함. 
sudo apt-get install libxml2-dev
sudo apt-get install mpich
transc c++버전을 설치 하기 전 neoParSA가 필요합니다
neoParSA는 cmake로 컴파일을 합니다. 

cmake 설치
sudo apt-get install cmake 

browse source에 neoParSA-1 (CMakeLink.txt가 있는 최상위 디렉토리)
browse build에 neoParSA-1/build (최상위 디렉토리 밑의 bulild 또는 bin 디렉토리)
File-delete cache 
configure
compiller : mpic++
아무 에러 없으면 generate

neoParSA-1/build 에 들어가서 make를 하면 
/usr/bin/ld: cannot find -lgsl
/usr/bin/ld: cannot find -lgslcblas
다음과 같은 에러메세지가 뜨는 데 
sudo apt-get install libgsl 또는 sudo apt-get install libgsl0ldbl
를 설치 해주면 해결


transc c++ 버전은 gcc와 g++로 compile을 합니다. 
하지만 libboost c++를 사용하기 때문에 따로 설치가 필요합니다. 
$ sudo apt-get install libboost-all-dev
or 
$ sudo aptitude install libboost-all-dev

의존성 문제는 해당 패키지를 설치해주면 됩니다. 

rcpp
$ sudo apt-get intall r-cran-cpp
or
in R 
> install.packages(“Rcpp”)

Makefile의 경로를 local에 맞게 바꿔줘야 합니다. 

MATLAB_DIR ?=/usr/local/matlab
RCPP_DIR   ?=/usr/lib/R/site-library/Rcpp/include
R_DIR      ?=/usr/bin/R
BOOST_DIR  ?=/usr/include/boost
PARSA_ROOT ?=/home/kang/src/neoParSA-1
PARSA_DIR = $(PARSA_ROOT)/parsa

LDLIBS += -L$(PARSA_ROOT)/build/lib -lparsa
LIBPARSA = $(PARSA_ROOT)/build/lib/libparsa.a

Rapi.R 파일을 조금 수정해 줘야 합니다. 
그냥 돌리면 오류가 남/ 오류나는 unfold의 path를 파일 수정을 통해 바꿔줘야 합니다. 

# 명령어 사용법
# scramble
transcpp인 xml 파일의 free parameter들을 랜덤하게 섞어주는 명령어
> scramble <원본 파일> <새로 섞인 파일>

원본파일의 내용 중 anneal="true" 로 되어있는 paramter의 value를 lim_low와 lim_high 사이의 랜덤한 값으로 설정해줌

# transcpp
free paramter들을 최적화 시켜주는 명령어 
scramble로 새로 생성한 파일을 사용한다. 
원본파일의 내용 중 anneal="true" 로 되어있는 paramter의 value를 lim_low와 lim_high 사이의 값 중 최적값을 찾아냄 
> transcpp <새로 섞인 파일>

학습과 관련 된 주요 parameter들은 아래와 같다

\<Root\>
 \<annealer_input init_T="100000" lambda="0.0001" init_loop="100000"/\>
 \<move interval="100" gain="3"/\>
 \<count_criterion freeze_crit="10" freeze_cnt="5"/\>
 \<mix adaptcoef="10"/\>
 \<lam tau="100" memLength_mean=".200" memLength_sd="100" criterion="10" freeze_cnt="5"/\>
 \<Mode\>
   \<ScoreFunction value="sse"/\> // sse, chisq, pdiff, rms, arkim, sss 가 가능함
   \<Competition value="false" window="500" shift="50"/\>
   \<NumThreads value="16"/\>
   \<SelfCompetition value="true"/\>
 
# unfold
