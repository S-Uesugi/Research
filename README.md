# Research
[実行コマンド]  
 nohup bash -c 'cat experiment.txt | xargs -P 5 -I {} bash exec.sh {}' & 
[experiment.txt]  
./experiment.out {実験する代数体のパス} {Kannanの埋め込み定数} {rotation の回数} {ブロックサイズ}
