# Research
[実行コマンド]  
 nohup bash -c 'cat experiment.txt | xargs -P 5 -I {} bash exec.sh {}' &   
[runfile.txt]  
./experiment.out {実験する代数体のパス} {Kannanの埋め込み定数} {rotation の回数} {ブロックサイズ}  
※ 実験する代数体のパスについて  
0 : exp1_131_1021  
1 : exp1_256_1021  
2 : original64_1021  
3 : original67_1021  
4 : original81_1021
