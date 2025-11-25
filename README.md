# Sigle-cell-RNA-seq
-- Lý thuyết về sigle cell RNA-seq có thể tham khảo:https://www.youtube.com/watch?v=jwSPTgF9ESQ&t=7087s và https://rnaseqcoban.github.io/R_tutorial/intro/
1.Các bước cơ bản trong thí nghiệm scRNA-seq.
- Phân lập tế bào từ mô.
- Phân tách RNA từ từng tế bào.
- Tổng hợp cDNA, nhân lên cDNA bằng phản ứng PCR hoặc nhân lên cDNA thông qua phản ứng phiên mã trong ống nghiệm (in vitro transcription) kết hợp với phiên mã ngược.
- Tạo thư viện DNA giải trình tự và giải trình tự.
- Phân tích kết quả.
<img width="814" height="579" alt="image" src="https://github.com/user-attachments/assets/f8e9ae3d-6ccc-4654-b0b3-cb03c8467062" />
2.Các bước cơ bản trong phân tích dữ liệu scRNA-seq.
Phân tích dữ liệu scRNAseq thường bắt đầu với “Expression matrix”. Trong Expression matrix, mỗi hàng biểu diễn một gene và mỗi cột biểu diễn một tế bào. Như vậy, mỗi ô biểu diễn mức độ biểu hiện của một gene trong một tế bào. Các bước xây dựng một “Expression matrix” bao gồm:
- Read Quality control: Đọc và kiểm tra chất lượng đoạn giải trình tự.
- Mapping Quality control: Kiểm tra độ chính xác của bước liên kết đoạn giải trình tự với reference genome dựa trên các chỉ số như: tỉ lệ giữa ribosome RNA và transposon RNA, tỉ lệ các đoạn nối đặc hiệu, độ sâu (read depth), đoạn đọc trải dài qua các điểm nối (splice junctions).
- Reads Quantification: Tính toán mức độ biểu hiện của mỗi gene trong một tế bào. Với phương pháp scRNA-seq sử dụng kỹ thuật “tag-base”, UMI có thể được dùng để tính số lượng tuyệt đối của một phân tử RNA.
  Công cụ xử lý đang tìm hiểu: Scanpy trong python 
3. Thách thức:
  - Khuếch đại không đồng đều cDNA (RT-PCR..)
  - Gen biểu hiện thấp ở tế bào và không xuất hiện ở các tế bào khác.
4. Phương pháp cải thiện:
Có rất nhiều phương pháp thí nghiệm, tuy nhiên, các phương pháp này có thể được phân nhóm dựa trên hai đặc điểm sau: (1) Cách định lượng RNA (2) Cách “bắt giữ” từng tế bào.

Có hai cách định lượng chính. Một là giải trình tự toàn bộ đoạn RNA (full-length). Hai là dùng “đầu dò” gắn đầu 3’ hoặc 5’, và chỉ giải trình tự đoạn gắn “đầu dò” (tag-based). Nên sử dụng phương pháp nào phụ thuộc vào mục đích tạo dữ liệu, vì cả hai đều có ưu và nhược điểm riêng. Về lí thuyết, “full-length” sẽ bao phổ (coverage) đồng đều đoạn phiên mã RNA. Tuy nhiên, thực tế cho thấy sự bao phổ (coverage) là không đồng đều. Nghĩa là có phần của đoạn phiên mã được bao phổ nhiều hơn so với các phần khác. Ưu điểm lớn nhất của phương pháp “tag-based” là có thêm trình tự nhận dạng phân tử đặc hiệu (Unique Molecular Identifier - UMI) trong đoạn mồi. Trình tự này là đặc trưng cho từng loại RNA khác nhau, vậy nên có thể được sử dụng trong định lượng. Tuy nhiên, nhược điểm của phương pháp này là bị giới hạn ở một đầu của RNA, và gây khó khăn trong việc phân biệt các isoform.(##có thể hiểu là gioongd đầu dò tagman trong RT-PCR về phương thức, nhược điểm là bao phủ không đồng đều nếu full-length và chỉ giớ hạn ở 1 đầu)

Khi nghĩ về phương pháp “bắt giữ” từng tế bào, chúng ta có thể xem xét ba sự lựa chọn sau: microwell, microfluidic, và droplet:

Microwell: sử dụng các phương pháp phân tách như FACS, pippette, laser capturing, để đưa từng tế bào vào trong một giếng có kích thước micro. Ưu điểm của microwell là có thể chụp hình ảnh của tế bào trong một giếng. Nhược điểm là quy mô nhỏ và khối lượng công việc lớn.

Microfulidic (Fluidigm C1): sử dụng nguyên lý thuỷ động lực học để dẫn dắt từng tế bào vào một buồng phản ứng. Ưu điểm là tính tự động hoá và quy mô cao hơn so với microwell. Tuy nhiên, hiệu suất thấp, chỉ có 10% số lượng tế bào được giữ lại ở các buồng phản ứng; và giá thành cho một đĩa “chip” quá đắt đỏ.

Droplet (10X Genomics, in-Drop): ý tưởng của phương pháp droplet là bao bọc một tế bào trong một giọt dầu đã có sẵn bead gắn đoạn mồi và nguyên liệu cho phản ứng tạo cDNA. Đây là phương pháp có hiệu suất và quy mô cao nhất trong ba loại.  (Giong phương pháp Droplet trong RT-PCR dùng 1 giọt dầu có sãn mồi để khuếch đại, trường hợp này là phiên mã ngược và khuếch đại CDNA trong cùng 1 phản ứng  ?? Trong phần luyện tập có bước loại bỏ double droplet)




