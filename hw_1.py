# 0. 라이브러리 로딩
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# 1. 전역 변수 설정
READ_COUNTS = 15
MIN_SEQ_LEN, MAX_SEQ_LEN = 201, 300      # 길이 범위
BASES = ("A", "G", "C", "T")     # 뉴클레오타이드
GLOBAL_SEED = 123                # 전체 재현성 시드

# 2. 개별 염기서열 생성 함수
def generate_random_sequence(seed, min_seq_len=MIN_SEQ_LEN, max_seq_len=MAX_SEQ_LEN, bases=BASES):  
    r = random.Random(seed)
    seq_len = r.randint(min_seq_len, max_seq_len)
    return "".join(r.choice(bases) for i in range(seq_len))

# 3. seed 값에 따른 염기서열 리스트 생성 함수
def generate_random_dna_sequences(global_seed, read_counts):
    rng = random.Random(global_seed)
    SEEDS = rng.sample(range(1, read_counts+1), read_counts)
    sequences = [generate_random_sequence(seed) for seed in SEEDS]
    return sequences
LIST_OF_SEQUENCES = generate_random_dna_sequences(GLOBAL_SEED, READ_COUNTS)

# 4. 염기서열 리스트 데이터프레임 생성
df_original = pd.DataFrame({  #수정
    'sequence_id': [f"Sequence_{i+1}" for i in range(len(LIST_OF_SEQUENCES))],  
    'sequence': LIST_OF_SEQUENCES,  
    'length': [len(seq) for seq in LIST_OF_SEQUENCES] 
})

# 5. 염기서열 리스트 txt 파일 생성
def save_txt(df, filename="sequences.txt"): 
    with open(filename, 'w') as file:  
            for i,row in df.iterrows():  
                file.write(f">{row['sequence_id']}\n") 
                file.write(f"{row['sequence']}\n") 
save_txt(df_original)

# 6. txt 파일 읽어오기
def read_txt(filename="sequences.txt"): 
    if not os.path.exists(filename): 
        print(f"파일 없음: {filename}") 
        return None 
    
    try: 
        sequences = []  
        headers = []  
        
        with open(filename, 'r') as file: 
            current_header = None 
            current_sequence = "" 
            
            for line in file: 
                line = line.strip()  
                if line.startswith('>'):
                    if current_header is not None:  
                        headers.append(current_header) 
                        sequences.append(current_sequence) 
                    current_header = line[1:]
                    current_sequence = ""  
                else: 
                    current_sequence += line
            
            if current_header is not None: 
                headers.append(current_header)
                sequences.append(current_sequence) 
        
        df_txt = pd.DataFrame({'sequence_id': headers, 'sequence': sequences}) 
        print(f"txt 파일 읽기 성공: {len(df_txt)}개 서열")
        return df_txt 
        
    except Exception as e: 
        print(f"txt 읽기 실패: {e}") 
        return None 

df_load=read_txt()

# 7. GC content percentage 계산 함수
def gc_function(sequence):
    gc_content = 100.0 * ( sequence.count('G') + sequence.count('C') ) / len(sequence)
    return gc_content

# 8. txt 파일에서 gc ratio 계산 결과 출력
df_read = df_load.copy()
df_read['gc_ratio'] = df_read['sequence'].apply(gc_function) 
df_read['length'] = df_read['sequence'].str.len()
df_read['sequence_number'] = df_read['sequence_id'].str.extract(r'Sequence_(\d+)').astype(int)

for i, row in df_read.iterrows():
  print(f"{row['sequence_id']}: Length = {row['length']}, GC Content = {row['gc_ratio']:.2f}%")

# 9. 시각화 
plt.figure(figsize=(8, 6))
plt.bar(df_read["sequence_number"], df_read["gc_ratio"])
plt.title('GC Content of DNA Sequences')
plt.xlabel('Sequence Number')
plt.ylabel('GC Content (%)')
plt.xticks(df_read["sequence_number"])
plt.yticks(range(0, 101, 20))
plt.ylim(0, 100)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
