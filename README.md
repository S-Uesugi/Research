# ソースコードの扱い方と各フォルダについて
初めて攻撃用のソースコードを実行する場合は必ずこのファイルを読んでください。  
このファイルでは攻撃コードの実行方法と実行に必要なフォルダとファイルについて説明します。  

1. **各フォルダについて**  
   攻撃用のコードを実行する場合、以下のフォルダ構成を用意して下さい。
   - **data** : Ring-LWEサンプルが格納されているフォルダです。サンプル生成のコードで生成したファイルをここに置いてください  
   （生成したファイルを一つのフォルダにまとめると複数の実験ができるようになります）。
   - **src** : 攻撃用のソースコードが入っているフォルダです。このフォルダ内には複数の実験が一括でできるようにrunfiles.txt,exec.shが格納されています。
   - **result** : 攻撃結果として出力されるファイルを格納するフォルダです。
2. **ソースコードの扱い方**  
[実行コマンド]  
 ```
 nohup bash -c 'cat experiment.txt | xargs -P 5 -I {} bash exec.sh {}' &   
```
  
[runfile.txt]  
./experiment.out {実験する代数体のパス} {Kannanの埋め込み定数} {rotation の回数} {ブロックサイズ}  
* 実験する代数体のパスについて  
0 : exp1_131_1021  
1 : exp1_256_1021  
2 : original64_1021  
3 : original67_1021  
4 : original81_1021
