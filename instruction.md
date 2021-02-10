# 惑大git入門編
# はじめに
この記事は基本的に[この公式ページ](https://git-scm.com/book/en/v2)で公開されている本をベースに書きました。もっと詳しく知りたい方は本を読んでください。
特にブランチの使い方と概念について理解するとgitと友達になれる気がしました。

# gitの概念

#### ・gitのデータの考え方

他のバージョン管理ツール -> ファイル毎に差分を記録していく
git -> **全てのファイルの状態をまとめてスナップショット的に記録する**
![Screen Shot 2020-06-22 at 22.21.52.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/c17812a8-fc07-ac68-1e9f-04daab9ca590.png)
[参照元](https://git-scm.com/book/en/v2/Getting-Started-What-is-Git%3F)

#### ・3つのステージ
<各ステージについて>
- working directory = 普段作業する場所
- staging area(index) = 次にコミットされるファイルが保存されてるgitディレクトリ
- .git directory(repository) = プロジェクトのメタデータがある場所(gitの最重要部)


![Screen Shot 2020-06-22 at 22.27.15.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/3db15bee-cb1b-dfd8-5827-b99125f3a116.png)


<作業の流れ>
1. 作業ディレクトリでファイルを編集
2. 編集したファイルのスナップショットをステージングエリアに追加

```
$ git add filename
```
3.コミットする(ステージングエリアにあるファイルを永久保持するスナップショットとしてgitディレクトリに保存する)

```
$ git commit -m "message"
```

# gitリポジトリ作成
主に2つ方法があります。
1. 既存のディレクトリでgitを初期化して作成

```
$ cd /好きな/ディレクトリ
$ git init
```

2. githubでリポジトリを作成してcloneする
githubアカウントを持ってない人は作ってください。
github上で新しいリポジトリを作成する。
clone or downloadをクリックしてそこにあるURLをコピー

```
$ git clone copied_URL
```

# git基本操作
## 基本サイクル

基本的にいつもやることは、

(0). gitの状態確認
(1). ファイル編集
(2). ステージングエリアに保存(git add)
(3). コミットする(git commit)

の繰り返しかと思います。
その後、どれかのブランチにpushという流れ(だと理解している)。

![Screen Shot 2020-06-22 at 23.22.56.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/d2ae47d9-8678-0361-012b-6642889208b4.png)

###### (0). gitの状態を確認します。

```
$ git status
```
すると,

- 自分が今どこのbranchにいるのか？
- Changes to be committed (次にコミットされる変更)
- Changes not staged for commit(まだステージングエリアにない変更)

などが分かります。
慣れていないうちは定期的にこのコマンドを使ってgitの状態を確認するといいと思います！

###### (1). ファイル編集
エディターで何か編集する。

###### (2). 編集したファイルをステージングエリアに！

```
$ git add filename
```

###### (3). ステージングエリアにある変更をコミットする

```
git commit -m "message"
```
addとcommitを一緒に行うコマンドもあります。

```
$ git commit -a -m 'some messages'
```

## 作業のやり直し！
誰でも間違うことはあります。

- 追加すべきファイルを追加せずコミットしてしまった or コミットメッセージが変

以下のようにすると、最終的なコミットは一つになる。

```
$ git commit -m 'initial commit'
$ git add forgotten_file
$ git commit --amend
```
- ステージしたファイルを取り消す

```
$ git reset HEAD unstaged_file
```

- ファイル変更の取り消し
これは個人的に一番使いそう。
色々変えてみたけど、全然うまくいかなかったので最初の状態に戻したい時です。

```
$ git checkout -- filename
```
これで直近のコミット状態に戻れます。

## 複数人でのリモート作業
gitは複数人で開発するのが普通ですよね。
リモートで作業する時の基本を簡単にまとめてみました。
イメージは以下のようになっています。

<img width="759" alt="Screen Shot 2020-06-23 at 21.12.12.png" src="https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/1889fb88-b08a-48e7-e3c8-f7004e7d74b5.png">

- どのリモートサーバーを設定したか確認

```
$ git remote -v
```

- fetchとmergeとpullの違い

fetch: リモートにある情報をローカルリポジトリに持ってくるだけ(作業ディレクトリのファイルは変わらない)
merge: ローカルリポジトリの情報を作業ディレクトリのファイルに反映
pull: fetch+mergeの作業を同時に行う

```
$ git fetch remote_name　#originなど
$ git pull origin branch_name
```

- リモートへpush

``git push remote_name branch_name``でできる。
**注意:**
**このpushは自分のローカル状態を最新にしてから、他に誰もリモートにpushしていない時しかできません。他の人がpushしていた場合には、一度pullしてローカル環境で調整してからでないとpushできない！**

## branchとは？
樹形図みたい感じでぼんやりと概念ではわかっていたつもりなのですが、自分は全く分かっていませんでした。
branch=ポインタです！！！
これを理解したらかなりスッキリしました。笑
まずブランチを理解するためには、gitがどうやってスナップショットを記録しているかを知る必要があります。
例えば、3つのファイルをaddしてcommitしたとする以下のような5つのオブジェクトがgitリポジトリに生成されます。コミットオブジェクトに編集者情報やツリーオブジェクト(真ん中の青いファイル)へのポインタ、そしてそのツリーオブジェクトにメタデータ(黄色の3つ)へのポインタが保存されています。
![Screen Shot 2020-06-23 at 22.25.58.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/f1cacecf-06d0-466a-db61-54c2849a8019.png)
[参照元](https://git-scm.com/book/ja/v2/Git-%E3%81%AE%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81%E6%A9%9F%E8%83%BD-%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81%E3%81%A8%E3%81%AF)

コミットを複数回行うと、新しいコミットに一つ前のコミットへのポインタが格納されます。
![Screen Shot 2020-06-23 at 22.34.29.png](https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/08074aa6-97e0-b36d-9295-e58509e3fcca.png)

**ここでブランチとは複数あるコミットの中でいずれかを指すポインタになります！**
それでは、最初masterしかない状態から新しいブランチを作成すると何が起こっているのかを図で見てみましょう。以下の作業を行います。
1. masterにコミット
2. testブランチを作成して移動
3. testブランチにコミット
4. 再びmasterに戻ってコミット

<img width="955" alt="Screen Shot 2020-06-23 at 23.05.34.png" src="https://qiita-image-store.s3.ap-northeast-1.amazonaws.com/0/199901/b341a9af-4f73-a200-46a6-6e88c55ab4c9.png">

最後に分岐したことが分かります。
また、ブランチがただポインタとして動いていることも理解できると思います。
HEADとは、自分が今どのブランチで作業しているかを示すポインタになります。つまり、gitはHEADというポインタを使ってどのブランチにコミットするか判断していることになりますね。
ここで注意です！！
gitでブランチを切り替えると、もちろん作業ディレクトリのファイルが変更されます。そのため、ブランチを切り替える前にファイルを変更したままにすると切り替えると怒られます。なので、**ブランチを切り替える時は、ちゃんと変更を今いるブランチにコミットしてから切り替えましょう！**

## ブランチの使い方(よくある例)
惑大でよく使うであろう作業内容をここに書きます。

#### ex.1) 誰かがorigin masterにpushした
1. まずgitの状態を確認！
```
$ git status
```
2. 次に自分のローカル環境を最新にするためにpullしましょう！
```
$ git checkout branch_name #まずはmasterに
$ git pull origin master #conflictしたらファイルを開いて解消する
```
3. コンフリクトしたらエディターで開いて解消する.-> [この記事へ](https://qiita.com/hkengo/items/f47b9f50ac2dca407d12)

#### ex.2) コードに変更を加える
1. まずは今回実装したい部分専用のbranchを作成します
```
$ git checkout feature/name #今回は自分のbranchから派生させる
$ git branch new_function　#ブランチ名はなんでも良い
$ git checkout new_function #new_branchに移動
```
2. コードに変更を加える->ステージング->コミット(ローカル)を繰り返す
```
$ # edit file
$ git add filename
$ git commit -m "message"
```
3. 上手く動くことが確認できたらfeature/nameにmergeする
```
$ git status #全て変更がコミットされてるか確認する。まだ変更したままだとエラー
$ git checkout feature/name #feature/nameブランチに戻る
$ git merge new_function #これでfeature/nameブランチが更新された
```
4. リモートのfeature/nameにpushする
```
$ git branch -a #feature/nameブランチにいるか確認。いなかったら移動する
$ git push origin feature/name
```
5. github上でpull requestする\
説明は[こちら](https://qiita.com/takamii228/items/80c0996a0b5fa39337bd)へ
6. 中村くんが承認してmasterにmergeする

#### ex.3) 前のあの状態に戻りたいと思ったら
- 指定したコミット時点の状態に戻る\
=作業ツリーを指定したコミット時点の状態にまで戻し、コミットを行う（コミットをなかったことにはせず、逆向きのコミットをすることで履歴を残す）
```
$ git log #これで戻りたいコミットのハッシュ値を確認してコピーする
$ git revert コミットのハッシュ値
```

- ファイル変更を取り消し(まだaddもcommitもしてない場合)
```
$ git checkout -- filename
```
補足: ファイル名だけでなく、ディレクトリや全部`.`も指定できる
- 直前のコミットをなかったことにする
```
$ git reset --hard HEAD^
```
`--hard`: コミット取り消した上でワークディレクトリの内容も書き換えたい場合に使用 \
`--soft`: ワークディレクトリの内容はそのままでコミットだけを取り消したい場合に使用 \
`HEAD^`: 直前のコミットを意味する\
`HEAD~{n}`: n個前のコミットを意味する(これの代わりにコミットのハッシュ値でもok)\

もっと詳しく知りたい方は以下のリンクへどうぞ。
- [ブランチとマージの基本<-これは必ず読むべき](https://git-scm.com/book/ja/v2/Git-%E3%81%AE%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81%E6%A9%9F%E8%83%BD-%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81%E3%81%A8%E3%83%9E%E3%83%BC%E3%82%B8%E3%81%AE%E5%9F%BA%E6%9C%AC)

- [リモートブランチのことをより深く理解できる<-ブランチの概念が分かっていればすんなり](https://git-scm.com/book/ja/v2/Git-%E3%81%AE%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81%E6%A9%9F%E8%83%BD-%E3%83%AA%E3%83%A2%E3%83%BC%E3%83%88%E3%83%96%E3%83%A9%E3%83%B3%E3%83%81)
- [共同開発の流れ<-さらっと読むだけ](https://git-scm.com/book/ja/v2/Git-%E3%81%A7%E3%81%AE%E5%88%86%E6%95%A3%E4%BD%9C%E6%A5%AD-%E5%88%86%E6%95%A3%E4%BD%9C%E6%A5%AD%E3%81%AE%E6%B5%81%E3%82%8C)


***

# git command集

### 基本操作

- 最新状態にアップデート(pull=fetch+merge)]

```
$ git pull origin master
```

- fetch: リモートの最新状態をローカルに持ってくる(ローカルのファイルには反映しない)

```
$ git fetch origin remote_branch_name
```
- merge: ローカルにあるリモートの最新状態(origin/master)をリモートのファイルに反映

```
$ git merge origin/master
```

- ファイル追加

```
$ git add filename
```

- コミットする

```
$ git commit -m “comment”
```
- ローカル->リモートのmasterにpush

```
$ git checkout master
$ git push origin master
```
- ローカルのあるブランチ->リモートのあるブランチにpush

```
$ git push origin <local branch name>:<remote branch name>
or
$ git push origin HEAD:<remote branch name>
```

- (ローカル)別ブランチのファイルを持ってくる

```
$ git chekout <brach name> <file name>
```

- (ローカル)別ブランチのディレクトリを持ってくる

```
$ git chekout <brach name> <directory name>
```
### 状態確認全般
- git状態の確認

```
$ git status
```

- 変更したけどまだステージしてない内容を確認

```
$ git diff
```

- ステージされた内容を確認

```
$ git diff --staged
```

- コミット履歴の確認

```
git log
```
オプションで
`-p`をつけると、反映された変更点を表示できる。これ便利。



### branch操作
- branch情報を見る

```
$ git branch -a
```

- branch 作成(新しくローカルで作る場合)

```
$ git branch new_branch_name
```
- branch 作成(リモートのブランチをローカルに持ってくる)


1. 新しいブランチをまずローカルに作成

```
$ git branch local_new_branch_name origin/branch_name
```
2. 新しいブランチに移動してpullする

```
$ git checkout local_new_branch_name
$ git pull origin branch_name
```

- branchの紐付けを確認する

```
$ git branch -vv
```

- ローカルのbranchを削除

```
$ git branch -D <branch name>
```

### 訂正・削除
- ローカルのファイル削除(正確には、ステージングエリアから削除し、コミットする+作業ディレクトリからも削除する)

```
$ git rm filename
```

- gitの中でファイル名変更

```
$ git mv file_from file_to
```

- リモートの状態に強制的に戻したい時

```
$ git fetch origin
$ git reset --hard origin/master
```

- リモートで削除されたブランチをローカルに反映させる

```
$ git remote prune origin
```

- 直前のコミットをなかったことにする

```
$ git reset --hard HEAD^
```
`--hard`: コミット取り消した上でワークディレクトリの内容も書き換えたい場合に使用 \
`--soft`: ワークディレクトリの内容はそのままでコミットだけを取り消したい場合に使用 \
`HEAD^`: 直前のコミットを意味する\
`HEAD~{n}`: n個前のコミットを意味する(これの代わりにコミットのハッシュ値でもok)\

- 指定したコミット時点の状態に戻る

=作業ツリーを指定したコミット時点の状態にまで戻し、コミットを行う（コミットをなかったことにはせず、逆向きのコミットをすることで履歴を残す）
```
$ git revert コミットのハッシュ値
```

- 直前のコミットの上書き

```
$ git commit --amend
```
