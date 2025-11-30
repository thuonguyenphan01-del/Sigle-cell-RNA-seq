# Sigle-cell-RNA-seq
File thực hành: Python (Scanpy và SCVI), R (Seurat).
Video hướng dẫn thực hành với R -studio: https://www.youtube.com/watch?v=5HBzgsz8qyk&t=463s.
script to perform standard workflow steps to analyze single cell RNA-Seq data.
data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1.
data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0.      
I. Lý thuyết về single cell RNA-seq
-- Lý thuyết về sigle cell RNA-seq có thể tham khảo:https://www.youtube.com/watch?v=jwSPTgF9ESQ&t=7087s và https://rnaseqcoban.github.io/R_tutorial/intro/.

1.Các bước cơ bản trong thí nghiệm scRNA-seq.
- Phân lập tế bào từ mô.
- Phân tách RNA từ từng tế bào.
- Tổng hợp cDNA, nhân lên cDNA bằng phản ứng PCR hoặc nhân lên cDNA thông qua phản ứng phiên mã trong ống nghiệm (in vitro transcription) kết hợp với phiên mã ngược.
- Tạo thư viện DNA giải trình tự và giải trình tự.
- Phân tích kết quả.
<img width="814" height="579" alt="image" src="https://github.com/user-attachments/assets/f8e9ae3d-6ccc-4654-b0b3-cb03c8467062" />.

2.Các bước cơ bản trong phân tích dữ liệu scRNA-seq.
Phân tích dữ liệu scRNAseq thường bắt đầu với “Expression matrix”. Trong Expression matrix, mỗi hàng biểu diễn một gene và mỗi cột biểu diễn một tế bào. Như vậy, mỗi ô biểu diễn mức độ biểu hiện của một gene trong một tế bào. Các bước xây dựng một “Expression matrix” bao gồm:
- Read Quality control: Đọc và kiểm tra chất lượng đoạn giải trình tự.
- Mapping Quality control: Kiểm tra độ chính xác của bước liên kết đoạn giải trình tự với reference genome dựa trên các chỉ số như: tỉ lệ giữa ribosome RNA và transposon RNA, tỉ lệ các đoạn nối đặc hiệu, độ sâu (read depth), đoạn đọc trải dài qua các điểm nối (splice junctions).
- Reads Quantification: Tính toán mức độ biểu hiện của mỗi gene trong một tế bào. Với phương pháp scRNA-seq sử dụng kỹ thuật “tag-base”, UMI có thể được dùng để tính số lượng tuyệt đối của một phân tử RNA.
  Công cụ xử lý đang tìm hiểu: Scanpy trong python.
  
3. Thách thức:
  - Khuếch đại không đồng đều cDNA (RT-PCR..).
  - Gen biểu hiện thấp ở tế bào và không xuất hiện ở các tế bào khác.
  - 
4. Phương pháp cải thiện:
- Có rất nhiều phương pháp thí nghiệm, tuy nhiên, các phương pháp này có thể được phân nhóm dựa trên hai đặc điểm sau: (1) Cách định lượng RNA (2) Cách “bắt giữ” từng tế bào.

- Có hai cách định lượng chính. Một là giải trình tự toàn bộ đoạn RNA (full-length). Hai là dùng “đầu dò” gắn đầu 3’ hoặc 5’, và chỉ giải trình tự đoạn gắn “đầu dò” (tag-based). Nên sử dụng phương pháp nào phụ thuộc vào mục đích tạo dữ liệu, vì cả hai đều có ưu và nhược điểm riêng. Về lí thuyết, “full-length” sẽ bao phổ (coverage) đồng đều đoạn phiên mã RNA. Tuy nhiên, thực tế cho thấy sự bao phổ (coverage) là không đồng đều. Nghĩa là có phần của đoạn phiên mã được bao phổ nhiều hơn so với các phần khác. Ưu điểm lớn nhất của phương pháp “tag-based” là có thêm trình tự nhận dạng phân tử đặc hiệu (Unique Molecular Identifier - UMI) trong đoạn mồi. Trình tự này là đặc trưng cho từng loại RNA khác nhau, vậy nên có thể được sử dụng trong định lượng. Tuy nhiên, nhược điểm của phương pháp này là bị giới hạn ở một đầu của RNA, và gây khó khăn trong việc phân biệt các isoform.(##có thể hiểu là gioongd đầu dò tagman trong RT-PCR về phương thức, nhược điểm là bao phủ không đồng đều nếu full-length và chỉ giớ hạn ở 1 đầu)

- Khi nghĩ về phương pháp “bắt giữ” từng tế bào, chúng ta có thể xem xét ba sự lựa chọn sau: microwell, microfluidic, và droplet:

- Microwell: sử dụng các phương pháp phân tách như FACS, pippette, laser capturing, để đưa từng tế bào vào trong một giếng có kích thước micro. Ưu điểm của microwell là có thể chụp hình ảnh của tế bào trong một giếng. Nhược điểm là quy mô nhỏ và khối lượng công việc lớn.

- Microfulidic (Fluidigm C1): sử dụng nguyên lý thuỷ động lực học để dẫn dắt từng tế bào vào một buồng phản ứng. Ưu điểm là tính tự động hoá và quy mô cao hơn so với microwell. Tuy nhiên, hiệu suất thấp, chỉ có 10% số lượng tế bào được giữ lại ở các buồng phản ứng; và giá thành cho một đĩa “chip” quá đắt đỏ.

- Droplet (10X Genomics, in-Drop): ý tưởng của phương pháp droplet là bao bọc một tế bào trong một giọt dầu đã có sẵn bead gắn đoạn mồi và nguyên liệu cho phản ứng tạo cDNA. Đây là phương pháp có hiệu suất và quy mô cao nhất trong ba loại.  (Giong phương pháp Droplet trong RT-PCR dùng 1 giọt dầu có sãn mồi để khuếch đại, trường hợp này là phiên mã ngược và khuếch đại CDNA trong cùng 1 phản ứng  ?? Trong phần luyện tập có bước loại bỏ double droplet).

II.Các bước phân tích dữ liệu:. 
1.Kiểm soát chất lượng dữ liệu (Cell Quality control).
- Một khi đã có dữ liệu, chúng ta cần kiểm tra và loại bỏ các tế bào có chất lượng dữ liệu kém. Sai sót trong việc loại bỏ các tế bào có dữ liệu kém có thể làm nhiễu dữ liệu và ảnh hưởng tới việc tìm ra các thông tin mang ý nghĩa sinh học.
- 
2.Chuẩn hoá dữ liệu (Normalizing data).
  - Bước chuẩn hóa dữ liệu luôn luôn cần việc xây dưng model. Tương tự nhue StandardScaler hay MinMax scaler .. thường hay dùng. Trong sigle cell RBA-seq thì thường dùng LogNormalization.
  - Ở phần này, chúng ta sẽ sử dụng phương pháp chuẩn hóa theo global-scaling LogNormalize để chuẩn hóa các phép đo gene expression cho mỗi tế bào với tổng lượng gene expression, nhân giá trị này với hệ số tỷ lệ (mặc định là 10.000) và biến đổi hàm log cho kết quả. Có rất nhiều phương pháp để chuẩn hóa dữ liệu, nhưng đây là phương pháp đơn giản và trực quan nhất. Phép chia cho tổng gene expression được thực hiện để thay đổi tất cả các số lượng gene expression thành một số đo tương đối. Kinh nghiệm cho thấy rằng các yếu tố kỹ thuật (ví dụ tỷ lệ bắt giữ được RNA, hiệu quả của phiên mã ngược) chịu trách nhiệm phần lớn cho sự thay đổi số lượng phân tử RNA trên mỗi tế bào.
  -  FindVariableFeatures là hàm tính gene expression trung bình và độ phân tán cho mỗi gen, đặt các gen này vào các bins, sau đó tính z-score cho sự phân tán trong mỗi bin. Điều này giúp kiểm soát mối quan hệ giữa độ biến thiên và gene expression trung bình.(Hiểu đơn giản là các gen có độ tương quan cao nhất với biểu hiện).
  - 
3.Giảm chiều dữ liệu (Dimensionality reduction). 
- PCA(đã học): Phương pháp giảm chiều dữ liệu phổ biến trong machine learning, dựa trên mô hình tuyến tính.Phương pháp này dựa trên quan sát rằng dữ liệu thường không phân bố ngẫu nhiên trong không gian mà thường phân bố gần các đường/mặt đặc biệt nào đó. PCA xem xét một trường hợp đặc biệt khi các mặt đặc biệt đó có dạng tuyến tính là các không gian con (subspace).Chỉ có thể locj các dữ liệu tuyến tính và bảo toàn phương sai cục bộ nhưng nó thường không nắm bắt được các mối quan hệ phi tuyến tính phức tạp. Vì vậy sẽ ưu tiên xử lý dữ kiệu thô trước rồi sẽ dùng 1 thuật toán giảm chiểu khác để xử lý dữ liệu phi tuyến tính.
- t-SNE(đã họ): Phù hợp để xử lý cá dự liệu phi tuyến nói cách khác là dữ liệu phân cụm. Giup trực quang hóa dữ liệu 2D và 3D.t-SNE emphasizes local relationships and clusters.
- UMAP: tương tự t-sne nhưng  t-SNE giữ lại cấu trúc cục bộ (có nghĩa là bạn không thể so sánh khoảng cách giữa các cụm) trong khi UMAP được cho là giữ lại cả cấu trúc cục bộ và toàn cục.UMAP focuses on complex, non-linear relationships.
  - Umap thường được xem là phiên bản cao cấp hơn của t-sne tuy nhiên nếu như mục tiêu của model là phân cụm riêng biêt và quan tâm đến mối quan hệ gần mà không cần quan tâm đến cấu trúc toạn cục thì nên dung t-sne.
4. Phân nhóm tế bào (Clustering).
  - Các thuật toán k-n cluster , Louran...
5.Phân tích khác biệt biểu hiện gen (Marker..)
6.Differential expression analysis.
Sau khi phân tế bào thành các nhóm khác nhau, câu hỏi đặt ra là điều gì đã làm cho mỗi nhóm khác biệt với các nhóm còn lại trong tập dữ liệu, và làm sao để chú thích chính xác các nhóm theo phân loại tế bào (tế bào máu, T-cells, B-cells,…) của chúng.

Có nhiều cách để thực hiện bước này:

Tìm kiếm các marker gene (gene có tính đại diện, nổi bật) làm tăng/giảm điều hòa (upregulation/downregulation) ở các loại tế bào đang quan tâm (so với phần còn lại của bộ dữ liệu).
So sách toàn bộ biểu hiện gene giữa các nhóm tế bào với nhau.
Sử dụng các phương pháp tự động so sánh các tế bào trên cơ sở dữ liệu về biểu hiện gene của các loại tế bào đã được công bố giúp dễ dàng hơn trong việc chú giải loại tế bào.
Các phương pháp tự động là một tiến bộ đầy hứa hẹn, nhưng vẫn chưa thể thay đổi được sự cẩn thận của con người.




