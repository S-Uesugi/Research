# Research
[実行コマンド] 
 nohup bash -c 'cat test.txt | xargs -P 5 -I {} bash exec.sh {}' &
