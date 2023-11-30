# Research
[実行コマンド] 
 nohup bash -c 'cat experiment.txt | xargs -P 5 -I {} bash exec.sh {}' &
